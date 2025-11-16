https://code.earthengine.google.com/001f24bad3c347c7a5a8f27d96bfc827

/* 

C. R. Baldwin - Dept. of Earth Sciences, University of Oxford
Mackenzie River surface area and Vitrekwa timeseries
For: Dasari et al. "Climate sensitive methane release from Arctic rivers in sediment-laden channels", November 2025 

Site coordinates (lat, long):
Vitrekwa-Peel tributary: 67.1665, -135.0688
Middle Channel start:	67.5785,	-134.0116
Middle Channel end: 68.4638,	-134.1263
East Channel start:	67.7887,	-134.1848
East Channel end: 68.3073,	-133.7333

*/


// 1. River surface area

// a. Define area of interest
var aoi = ee.Geometry.Polygon(
  [[[-134.5911, 68.4828],
    [-134.5911, 67.5168],
    [-133.6545, 67.5168],
    [-133.6545, 68.4828]]], null, false);

Map.addLayer(ee.Image().paint(aoi, 0, 1), {palette: 'black'}, 'AOI outline');
Map.centerObject(aoi, 8);

// b. user parameters
var RSAstart = 2020;
var RSAend   = 2025; // this is for a 5-yr averaged composite
var scale = 10; // m/px
var sentinel_collection = 'COPERNICUS/S2_SR_HARMONIZED'; // specific satellite image collection we want to use
var sampledates = ['2023-09-17', '2024-06-21']; // this is for specific sampling trip dates as well as the composite

// c. functions
function maskS2Clouds(image) {
  var qa = image.select('QA60');
  var mask = qa.bitwiseAnd(1 << 10).eq(0)
               .and(qa.bitwiseAnd(1 << 11).eq(0));
  return image.updateMask(mask);
}

function addNDWI(image) {
  var ndwi = image.normalizedDifference(['B3', 'B8']).rename('NDWI');
  return image.addBands(ndwi);
}


// d. set up satellite image collection using above functions 
var satimgs = ee.ImageCollection(sentinel_collection)
  .filterBounds(aoi);

var years = ee.List.sequence(RSAstart, RSAend);

var summer_collection = ee.ImageCollection(
  years.iterate(function(y, acc) {
    var start = ee.Date.fromYMD(y,6,10); // earliest ice-free image across all years
    var end   = ee.Date.fromYMD(y,9,20); // latest ice-free image
    
    var yearly = satimgs
      .filterDate(start,end)
      .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 15))
      .map(maskS2Clouds)
      .map(addNDWI);
    
    return ee.ImageCollection(acc).merge(yearly);
  }, ee.ImageCollection([]))
);

// e. create composite images to get a 5-year average of channel area
var ndwi_composite = summer_collection.select('NDWI').median().clip(aoi);

Map.addLayer(ndwi_composite, {min:-1, max:1, palette:['white','lightblue','blue']}, 'NDWI Composite');

Export.image.toDrive({
  image: ndwi_composite.toFloat(),
  description: 'NDWI_Composite_2020_2025',
  folder: 'GEE',
  fileNamePrefix: 'NDWI_Composite_2020_2025',
  region: aoi,
  scale: scale,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});


var rgb_composite = summer_collection
  .select(['B2','B3','B4']) // only reflectance bands, exclude NDWI
  .median()
  .clip(aoi)
  .toFloat(); // ensure consistent type

// Map.addLayer(rgb_composite, {bands:['B4','B3','B2'], min:0, max:3000}, 'RGB Composite');

Export.image.toDrive({
  image: rgb_composite,
  description: 'RGB_Composite_2020_2025',
  folder: 'GEE',
  fileNamePrefix: 'RGB_Composite_2020_2025',
  region: aoi,
  scale: scale,
  crs: 'EPSG:4326',
  maxPixels: 1e13
});

// f. date-specific images to align with sampling trips
sampledates.forEach(function(dateStr) {
  var start = ee.Date(dateStr);
  var end   = start.advance(1, 'day');

  // NDWI
  var ndwi_img = satimgs
    .filterDate(start, end)
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 15))
    .map(maskS2Clouds)
    .map(addNDWI)
    .select('NDWI')
    .mosaic()
    .clip(aoi);

  if (ndwi_img.bandNames().size().getInfo() === 0) {
    print('No NDWI data for', dateStr);
    return;
  }

  // Map.addLayer(ndwi_img, {min:-1, max:1, palette:['white','lightblue','blue']}, 'NDWI ' + dateStr);

  Export.image.toDrive({
    image: ndwi_img.toFloat(),
    description: 'NDWI_' + dateStr,
    folder: 'GEE',
    fileNamePrefix: 'NDWI_' + dateStr,
    region: aoi,
    scale: scale,
    crs: 'EPSG:4326',
    maxPixels: 1e13
  });

  // RGB
  var rgb_img = satimgs
    .filterDate(start, end)
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 15))
    .map(maskS2Clouds)
    .median()
    .select(['B2','B3','B4']) // ensure only reflectance bands
    .clip(aoi)
    .toFloat();

  if (rgb_img.bandNames().size().getInfo() === 0) {
    print('No RGB data for', dateStr);
    return;
  }

  // Map.addLayer(rgb_img, {bands:['B4','B3','B2'], min:0, max:3000}, 'RGB ' + dateStr);

  Export.image.toDrive({
    image: rgb_img,
    description: 'RGB_' + dateStr,
    folder: 'GEE',
    fileNamePrefix: 'RGB_' + dateStr,
    region: aoi,
    scale: scale,
    crs: 'EPSG:4326',
    maxPixels: 1e13
  });
});


// end sequence, now run Tasks and load into ArcGIS for raster to polygon process

// 2. Vitrekwa tributary split. First define dates (scouted using Copernicus Browser)
// event happened around 2025-05-24/27, new channel fully open by 2025-06-15

// a. define dates and AOI
var VK_0617 = '2017-06-20'; // farthest back S2 will go for summer images
var VK_0618 = '2018-06-20';
var VK_0619 = '2019-06-18';
var VK_0620 = '2020-06-14';
var VK_0621 = '2021-06-17';
var VK_0622 = '2022-06-21';
var VK_0623 = '2023-06-17';
var VK_0624 = '2024-06-13';
var VK_0525 = '2025-05-27'; // ice has melted and channel has just breached bank
var VK_0625 = '2025-06-16'; //another June image for consistency
var VK_0725 = '2025-07-28'; // mid-summer, can see sediment accumuation clearly
var VK_0925 = '2025-09-04'; // latest ice-free image showing most recent channel morphology

var VKdates = [
  {name: 'VK_0617', date: VK_0617},
  {name: 'VK_0618', date: VK_0618},
  {name: 'VK_0619', date: VK_0619},
  {name: 'VK_0620', date: VK_0620},
  {name: 'VK_0621', date: VK_0621},
  {name: 'VK_0622', date: VK_0622},
  {name: 'VK_0623', date: VK_0623},
  {name: 'VK_0624', date: VK_0624},
  {name: 'VK_0525', date: VK_0525},
  {name: 'VK_0625', date: VK_0625},
  {name: 'VK_0725', date: VK_0725},
  {name: 'VK_0925', date: VK_0925}
];

var VKaoi = ee.Geometry.Polygon(
  [[[-135.08611381335444, 67.17522713188508],
    [-135.08611381335444, 67.15444213890186],
    [-134.9635475902099,  67.15444213890186],
    [-134.9635475902099,  67.17522713188508]]], null, false
);


function export_VK(entry) {
  var date = ee.Date(entry.date);

  var VKimg = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filterBounds(VKaoi)
    .filterDate(date, date.advance(1, 'day'))
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 15))
    .map(maskS2Clouds)
    .map(addNDWI)
    .median() //sometimes there are two or three image tiles for a specific date so take a median of all tiles
    .clip(VKaoi);


  // RGB 
  var rgb = VKimg.select(['B4','B3','B2']).toFloat();
  // Map.addLayer(rgb, {bands:['B4','B3','B2'], min:0, max:3000}, entry.name + ' RGB');

  Export.image.toDrive({
    image: rgb,
    description: entry.name + '_RGB',
    folder: 'GEE',
    fileNamePrefix: entry.name + '_RGB',
    region: VKaoi,
    scale: 10,
    crs: 'EPSG:4326',
    maxPixels: 1e13
  });

  // NDWI - left out for now
  // var ndwi = VKimg.select('NDWI').toFloat();
  // Map.addLayer(ndwi, {min:-1, max:1, palette:['white','lightblue','blue']}, entry.name + ' NDWI');

  // Export.image.toDrive({
  //  image: ndwi,
  //  description: entry.name + '_NDWI',
  //  folder: 'GEE',
  //  fileNamePrefix: entry.name + '_NDWI',
  //  region: VKaoi,
  //  scale: 10,
  //  crs: 'EPSG:4326',
  //  maxPixels: 1e13
  // });


}

// Run for all dates
VKdates.forEach(export_VK);
