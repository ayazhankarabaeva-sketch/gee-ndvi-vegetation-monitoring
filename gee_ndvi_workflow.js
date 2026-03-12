var zona = ee.FeatureCollection([
  ee.Feature(geometry, {id: 1}),
  ee.Feature(geometry2, {id: 2})
]);

var roi = zona.geometry();
Map.centerObject(roi, 13);
Map.addLayer(zona, {color: '000000'}, 'zona', true);

var tzStart = ee.Date('2025-05-19');
var tzEnd   = ee.Date('2025-09-18');

var targets = ee.List([
  '2025-05-19',
  '2025-06-19',
  '2025-07-19',
  '2025-08-19',
  '2025-09-18'
]);

var windowDays  = 10;
var maxCloudPct = 20;
var TOPK        = 12;
var PNG_DIM     = 2200;

var rgbVis = {
  bands: ['B4', 'B3', 'B2'],
  min: 0.01,
  max: 0.30,
  gamma: 1.25
};

var naturalPalette5 = [
  '#b7a58a',
  '#dce9a5',
  '#b7db7a',
  '#6fbf6a',
  '#2b8c4b'
];

var zoneVisNatural = {
  min: 1,
  max: 5,
  palette: naturalPalette5
};

var NDVI_BREAKS = [0.2, 0.4, 0.6, 0.8];
var EVI_BREAKS  = [0.1, 0.2, 0.35, 0.5];

var field1 = ee.Feature(zona.filter(ee.Filter.eq('id', 1)).first());
var field2 = ee.Feature(zona.filter(ee.Filter.eq('id', 2)).first());

function sclClearMask(img) {
  var scl = img.select('SCL');
  return scl.neq(3)
    .and(scl.neq(8))
    .and(scl.neq(9))
    .and(scl.neq(10))
    .and(scl.neq(11));
}

function maskS2sr(img) {
  var clear = sclClearMask(img);
  var dataMask = img.select('B8').mask();
  return img.updateMask(clear).updateMask(dataMask);
}

function addIndices(img) {
  var nir  = img.select('B8').multiply(0.0001);
  var red  = img.select('B4').multiply(0.0001);
  var blue = img.select('B2').multiply(0.0001);

  var ndvi = nir.subtract(red).divide(nir.add(red)).rename('NDVI');

  var evi = img.expression(
    '2.5 * (nir - red) / (nir + 6 * red - 7.5 * blue + 1)',
    {nir: nir, red: red, blue: blue}
  ).rename('EVI');

  return img.addBands([ndvi, evi]);
}

function fractionInGeom(mask01, geom, scaleVal) {
  var v = mask01.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: geom,
    scale: scaleVal,
    bestEffort: true,
    tileScale: 8,
    maxPixels: 1e13
  }).values().get(0);
  return ee.Number(v);
}

function scoreCandidate(img, targetDate) {
  img = ee.Image(img);

  var ts = ee.Date(img.get('system:time_start'));
  var diff = ts.difference(targetDate, 'day').abs();

  var dataMask = img.select('B8').mask();
  var clearMask = sclClearMask(img).and(dataMask);

  var coverFrac = fractionInGeom(dataMask, roi, 20);
  var clearFrac = fractionInGeom(clearMask, roi, 20);
  var cloudPct = ee.Number(img.get('CLOUDY_PIXEL_PERCENTAGE'));

  return img.set({
    diff_days: diff,
    cover_frac: coverFrac,
    clear_frac: clearFrac,
    roi_cloud_frac: ee.Number(1).subtract(clearFrac),
    cloud_pct: cloudPct
  });
}

function zonesByBreaks(x, breaks) {
  breaks = ee.List(breaks);
  var b1 = ee.Number(breaks.get(0));
  var b2 = ee.Number(breaks.get(1));
  var b3 = ee.Number(breaks.get(2));
  var b4 = ee.Number(breaks.get(3));

  var z = ee.Image(0)
    .where(x.lt(b1), 1)
    .where(x.gte(b1).and(x.lt(b2)), 2)
    .where(x.gte(b2).and(x.lt(b3)), 3)
    .where(x.gte(b3).and(x.lt(b4)), 4)
    .where(x.gte(b4), 5);

  return z.updateMask(x.mask());
}

function pickSceneForTarget(targetStr) {
  targetStr = ee.String(targetStr);
  var targetDate = ee.Date(targetStr);

  var start = targetDate.advance(-windowDays, 'day');
  var end   = targetDate.advance(windowDays, 'day');

  start = ee.Date(ee.Algorithms.If(start.millis().lt(tzStart.millis()), tzStart, start));
  end   = ee.Date(ee.Algorithms.If(end.millis().gt(tzEnd.millis()), tzEnd, end));

  var base = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filterBounds(roi)
    .filterDate(start, end.advance(1, 'day'))
    .filter(ee.Filter.lte('CLOUDY_PIXEL_PERCENTAGE', maxCloudPct));

  var prelim = base.sort('system:time_start');
  var topCol = ee.ImageCollection(prelim.toList(TOPK)).map(function(im) {
    return scoreCandidate(im, targetDate);
  });

  var n = topCol.size();

  var best = ee.Image(
    topCol
      .sort('cloud_pct')
      .sort('diff_days')
      .sort('clear_frac', false)
      .sort('cover_frac', false)
      .first()
  );

  var ts = best.get('system:time_start');
  var imageDate = ee.Date(ts).format('YYYY-MM-dd');
  var proj10 = best.select('B4').projection();
  var roiMask = ee.Image(1).clip(roi).mask();
  var bestMasked = maskS2sr(best);

  var rgb = bestMasked.select(['B4', 'B3', 'B2'])
    .multiply(0.0001)
    .clip(roi)
    .updateMask(roiMask);

  var idx = addIndices(bestMasked).select(['NDVI', 'EVI'])
    .setDefaultProjection(proj10)
    .reproject(proj10)
    .clip(roi)
    .updateMask(roiMask);

  var packed = rgb.addBands(idx).set({
    target_date: targetStr,
    image_date: imageDate,
    n_candidates: n,
    cover_frac: best.get('cover_frac'),
    clear_frac: best.get('clear_frac'),
    roi_cloud_frac: best.get('roi_cloud_frac'),
    cloud_pct: best.get('cloud_pct'),
    diff_days: best.get('diff_days')
  }).set('system:time_start', ts);

  var empty = ee.Image(0).rename('B4')
    .addBands(ee.Image(0).rename('B3'))
    .addBands(ee.Image(0).rename('B2'))
    .addBands(ee.Image(0).rename('NDVI'))
    .addBands(ee.Image(0).rename('EVI'))
    .updateMask(ee.Image(0))
    .set({
      target_date: targetStr,
      image_date: targetStr,
      n_candidates: 0,
      cover_frac: 0,
      clear_frac: 0,
      roi_cloud_frac: 1,
      cloud_pct: null,
      diff_days: null
    })
    .set('system:time_start', targetDate.millis());

  return ee.Image(ee.Algorithms.If(n.gt(0), packed, empty));
}

var scenes = ee.ImageCollection(targets.map(pickSceneForTarget));
var scenesList = scenes.toList(5);

var ticks = ee.List.sequence(0, 4).map(function(i) {
  return ee.Image(scenesList.get(i)).get('image_date');
});

print('Targets:', targets);
print('Chosen image dates:', ticks);
print('cover_frac:', scenes.aggregate_array('cover_frac'));
print('clear_frac:', scenes.aggregate_array('clear_frac'));
print('roi_cloud_frac:', scenes.aggregate_array('roi_cloud_frac'));
print('diff_days:', scenes.aggregate_array('diff_days'));
print('tile cloud%:', scenes.aggregate_array('cloud_pct'));

function makeLegendInner(title, labels, palette) {
  var panel = ui.Panel({style: {padding: '0px', margin: '0 0 10px 0'}});
  panel.add(ui.Label({value: title, style: {fontWeight: 'bold', margin: '0 0 6px 0'}}));

  for (var i = 0; i < labels.length; i++) {
    var box = ui.Label('', {
      backgroundColor: palette[i],
      padding: '8px',
      margin: '0 6px 4px 0'
    });
    var lab = ui.Label(labels[i], {
      margin: '0 0 4px 0',
      fontSize: '12px'
    });
    panel.add(ui.Panel([box, lab], ui.Panel.Layout.flow('horizontal')));
  }
  return panel;
}

var ndviLegendLabels = [
  'Zone 1: < 0.20 — soil / almost no vegetation',
  'Zone 2: 0.20–0.40 — sparse vegetation',
  'Zone 3: 0.40–0.60 — medium vegetation',
  'Zone 4: 0.60–0.80 — good vegetation',
  'Zone 5: > 0.80 — very dense vegetation'
];

var eviLegendLabels = [
  'Zone 1: < 0.10 — soil / almost no vegetation',
  'Zone 2: 0.10–0.20 — sparse vegetation',
  'Zone 3: 0.20–0.35 — medium vegetation',
  'Zone 4: 0.35–0.50 — good vegetation',
  'Zone 5: > 0.50 — very dense vegetation'
];

var legendsBox = ui.Panel({
  style: {
    position: 'bottom-right',
    padding: '8px',
    width: '380px'
  }
});

legendsBox.add(makeLegendInner('NDVI zones', ndviLegendLabels, naturalPalette5));
legendsBox.add(makeLegendInner('EVI zones', eviLegendLabels, naturalPalette5));
Map.add(legendsBox);

for (var i = 0; i < 5; i++) {
  var img = ee.Image(scenesList.get(i));

  var d = ee.String(img.get('image_date')).getInfo();
  var t = ee.String(img.get('target_date')).getInfo();
  var cf = ee.Number(img.get('clear_frac')).format('%.2f').getInfo();
  var cov = ee.Number(img.get('cover_frac')).format('%.2f').getInfo();
  var roiCloud = ee.Number(img.get('roi_cloud_frac')).format('%.2f').getInfo();

  var ndvi = img.select('NDVI');
  var evi  = img.select('EVI');

  var ndviZones = zonesByBreaks(ndvi, NDVI_BREAKS);
  var eviZones  = zonesByBreaks(evi, EVI_BREAKS);

  Map.addLayer(
    img.select(['B4', 'B3', 'B2']),
    rgbVis,
    'RGB | target=' + t + ' | image=' + d + ' | cover=' + cov + ' | clear=' + cf + ' | roi_cloud=' + roiCloud,
    false
  );

  Map.addLayer(
    ndviZones,
    zoneVisNatural,
    'NDVI zoned | ' + d,
    false
  );

  Map.addLayer(
    eviZones,
    zoneVisNatural,
    'EVI zoned | ' + d,
    false
  );
}

function makeChart(field, label) {
  var fc = ee.FeatureCollection(scenesList.map(function(img) {
    img = ee.Image(img);

    var stats = img.select(['NDVI', 'EVI']).reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: field.geometry(),
      scale: 10,
      bestEffort: true,
      tileScale: 4,
      maxPixels: 1e13
    });

    return ee.Feature(null, {
      date: img.get('image_date'),
      NDVI: stats.get('NDVI'),
      EVI: stats.get('EVI')
    });
  }));

  var chart = ui.Chart.feature.byFeature(fc, 'date', ['NDVI', 'EVI'])
    .setChartType('LineChart')
    .setOptions({
      title: label + ' — NDVI & EVI mean',
      hAxis: {
        title: 'Image date',
        slantedText: true,
        slantedTextAngle: 45
      },
      vAxis: {title: 'Index'},
      legend: {position: 'bottom'},
      lineWidth: 2,
      pointSize: 6,
      chartArea: {left: 70, right: 20, top: 40, bottom: 140}
    });

  print(chart);
  print('Table — ' + label + ':', fc);
}

makeChart(field1, 'Field 1');
makeChart(field2, 'Field 2');

function printPngLinks(img) {
  img = ee.Image(img);

  var t = ee.String(img.get('target_date')).getInfo();
  var d = ee.String(img.get('image_date')).getInfo();

  var ndvi = img.select('NDVI');
  var evi  = img.select('EVI');

  var ndviZones = zonesByBreaks(ndvi, NDVI_BREAKS);
  var eviZones  = zonesByBreaks(evi, EVI_BREAKS);

  var rgbPng  = img.select(['B4', 'B3', 'B2']).visualize(rgbVis).clip(roi);
  var ndviPng = ndviZones.visualize(zoneVisNatural).clip(roi);
  var eviPng  = eviZones.visualize(zoneVisNatural).clip(roi);

  var params = {
    region: roi,
    dimensions: PNG_DIM,
    format: 'png'
  };

  print('PNG RGB  | target=' + t + ' | image=' + d, rgbPng.getThumbURL(params));
  print('PNG NDVI | target=' + t + ' | image=' + d, ndviPng.getThumbURL(params));
  print('PNG EVI  | target=' + t + ' | image=' + d, eviPng.getThumbURL(params));
}

print('--- PNG DOWNLOAD LINKS ---');

for (var j = 0; j < 5; j++) {
  printPngLinks(ee.Image(scenesList.get(j)));
}
