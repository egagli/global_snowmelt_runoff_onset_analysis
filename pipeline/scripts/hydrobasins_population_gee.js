// Recovered recipe for data/geometries/Hydrobasins_L5_Population_Global.geojson
// (298 MB, exported April 2025). Found in Eric's Earth Engine account on
// 2026-09-02; saved here verbatim so the file is regenerable. It was run in the
// Earth Engine Code Editor (JavaScript API) and exported to Google Drive.
//
// What it does: sums the GPW v4.11 population-count raster (CIESIN, most recent
// epoch = 2020, native 30 arc-second grid, ~927 m) inside every HydroBASINS v1
// level-5 polygon (WWF/HydroSHEDS/v1/Basins/hybas_5 — the same HydroBASINS
// geometries and PFAF_ID / HYBAS_ID keys as the BasinATLAS lev05 layer the
// pipeline rasterizes) and keeps HYBAS_ID, PFAF_ID, total_population, geometry.
//
// Planned after the v10 campaign: port to a
// Python entry point (ee API via gsro_analysis.settings.initialize_earthengine)
// at HydroBASINS level 6 — the basin key stored in the partials since the
// 2026-09-02 decision — so population is available per level-6 basin and sums
// exactly to level 5 by PFAF_ID prefix. Until then this file is the recipe.

// Script to aggregate population estimates into HydroBASINS level 5 product
// with only the required columns: HYBAS_ID, PFAF_ID, total_population, and geometry

// Load the global population count dataset
var populationDataset = ee.ImageCollection('CIESIN/GPWv411/GPW_Population_Count');
// Get the most recent population data
var populationImage = populationDataset.sort('system:time_start', false).first();
var population = populationImage.select('population_count');

// Get the native resolution of the population data
// Always use native resolution for population calculations to maintain accuracy
var nativeProjection = population.projection();
var nativeScale = nativeProjection.nominalScale();
print('Native population data resolution (meters):', nativeScale);

// Load the HydroBASINS level 5 dataset
var hydrobasins = ee.FeatureCollection('WWF/HydroSHEDS/v1/Basins/hybas_5');
print('Number of basins:', hydrobasins.size());

// Function to calculate population sum for each basin and keep only the required columns
var calculatePopulation = function(basins) {
  return population.reduceRegions({
    collection: basins,
    reducer: ee.Reducer.sum(),
    scale: nativeScale,
    tileScale: 16  // Prevent computation timeout
  }).map(function(feature) {
    // Keep only the required properties: HYBAS_ID, PFAF_ID, total_population
    return ee.Feature(feature.geometry(), {
      'HYBAS_ID': feature.get('HYBAS_ID'),
      'PFAF_ID': feature.get('PFAF_ID'),
      'total_population': feature.get('sum')
    });
  });
};

// Calculate population for all basins
var basinsWithPopulation = calculatePopulation(hydrobasins);

// Display the results
Map.setCenter(0, 20, 3);

// Layer 1 (bottom): Original population raster
// Visualization parameters for population raster
var popVis = {
  min: 0,
  max: 1000,
  palette: [
    'ffffe7',
    '86a192',
    '509791',
    '307296',
    '2c4484',
    '000066'
  ]
};
// Add the population layer first (will be on the bottom)
Map.addLayer(population, popVis, 'Population Count', true);

// Layer 2 (middle): Population aggregated into basins
// Convert feature collection to image for visualization
var populationImage = basinsWithPopulation.reduceToImage({
  properties: ['total_population'],
  reducer: ee.Reducer.first()
});

// Viridis-like color palette (perceptually uniform, colorblind-friendly)
var viridisColors = ['440154', '414487', '2a788e', '22a884', '7ad151', 'fde725'];

// Add the aggregated population layer (will be in the middle)
Map.addLayer(populationImage, {
  min: 0,
  max: 1000000,
  palette: viridisColors
}, 'Hydrobasins Level 5 (Population)', true);

// Layer 3 (top): Basin outlines
// Add basin outlines last (will be on top)
Map.addLayer(hydrobasins.style({
  color: 'white',
  fillColor: '00000000', // Transparent fill
  width: 1
}), null, 'Hydrobasins Level 5 (Outlines)', true);

// Export the results for the entire world
Export.table.toDrive({
  collection: basinsWithPopulation,
  description: 'Hydrobasins_L5_Population_Global',
  fileFormat: 'GeoJSON'
});

// OPTIONAL: Regional Processing
// If the global computation times out, use these regional exports

// Process by continent using HYBAS_ID first digit (commented out by default)
/*
// Function to export by continent ID
var exportByContinent = function(continentId) {
  // Filter basins by the first digit of HYBAS_ID
  var continentBasins = hydrobasins.filter(
    ee.Filter.rangeContains('HYBAS_ID', continentId * 1000000000, (continentId + 1) * 1000000000 - 1));
  var continentBasinsWithPop = calculatePopulation(continentBasins);

  Export.table.toDrive({
    collection: continentBasinsWithPop,
    description: 'Hydrobasins_L5_Population_Continent_' + continentId,
    fileFormat: 'GeoJSON'
  });
};

// Export each continent separately
for(var i = 1; i <= 9; i++) {
  exportByContinent(i);
}
*/
