/*                                              
======RUSLE Model for the case of Crete==========
*/          

//Bounding Box of area in Crete, Greece
var box = ee.FeatureCollection("users/alexandridisvasileios/box_wgs84");

//====== C-factor=============//
//Mean Sentinel-2 summer imagery for the year 2019//
function maskS2clouds(image) {
  var qa = image.select('QA60');

  // Bits 10 and 11 are clouds and cirrus, respectively.
  var cloudBitMask = 1 << 10;
  var cirrusBitMask = 1 << 11;

  // Both flags should be set to zero, indicating clear conditions.
  var mask = qa.bitwiseAnd(cloudBitMask).eq(0)
      .and(qa.bitwiseAnd(cirrusBitMask).eq(0));

  return image.updateMask(mask).divide(10000);
}
var S2L2A = ee.ImageCollection('COPERNICUS/S2_SR')
                  .select(['B4','B8','QA60'])
                  .filterDate('2020-06-01','2020-08-31')
                  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE',1))
                  .map(maskS2clouds)
                  .map(function(image){return image.clip(box.geometry())}).median();
//C-factor with NDVI
var findndvi = function(image){return image.addBands(image.normalizedDifference(['B8', 'B4']).clip(box.geometry()).rename('ndvi'));}
var ndvi_img2020 = findndvi(S2L2A).select(['ndvi']);
var cfactor_ndvi = ((ndvi_img2020.multiply(-2)).divide(ndvi_img2020.multiply(-1).add(1))).exp();

//====== R-factor=============//
//-----ERA5 daily percipitation data-----//
var precip = function(image){return image.expression('P*1000', {P: image.select('total_precipitation')}).set('system:time_start',image.get('system:time_start')).float()};

//Positive values -- bigger than 0
var positive = function(image)
{ var mymask = image.gte(0); 
  var updated = image.updateMask(mymask);
  var newa = updated.unmask(0).clip(box.geometry()); 
  return newa
}  

//Function to calculate R-factor per month for EACH individual year
function calc(startDate,endDate,box){
    var mymonths = ee.List.sequence(1,12);
    
    var ERA5_daily = ee.ImageCollection("ECMWF/ERA5/DAILY")
                  .filterDate(startDate,endDate)
                  .select('total_precipitation')//select precipitation
                  .map(precip)//convert to mm
                  .map(function(image){return image.clip(box.geometry())}); //Clips data based on 'aoi';
                  
    //R--index -- Apply precip function, based on de Santos Loureiro
    //Bigger than 10, the pixel values of percipitation
    var bigger10 = function(image)
    {   var mymask = image.gte(10); 
      var updated = image.updateMask(mymask);
      var newa = updated.unmask(0).clip(box.geometry()); 
      return newa
    }
    //Count the values bigger than 10
    var countD10 = function(image)
    { var mymask = image.lt(10);  
      var updated = image.updateMask(mymask);
      var newa = updated.unmask(1).clip(box.geometry()); 
      return newa
    }
    
    //get ERA5 bigger than 10, get counter with ERA5 bigger than 10
    var ERA_bigger = ERA5_daily.map(bigger10);
    var ERA_counter = ERA_bigger.map(countD10);
    //sum them!
    var bigger_sum = ERA_bigger.sum();
    var counter_sum = ERA_counter.sum();
    
    // Group by month R10, and then reduce within groups by sum();
    var bymonthR10 = ee.ImageCollection.fromImages(
        mymonths.map(function(m){
        var mm = ERA5_daily.filter(ee.Filter.calendarRange(m,m,'month'))
                          .select('total_precipitation')
                          .map(bigger10).sum() // NOT applied in a client-side operation
                          .set('month',m).rename('R10');
        return mm;
    }));
  
    // Group by month D10, and then reduce within groups by sum();
    var bymonthD10 = ee.ImageCollection.fromImages(
        mymonths.map(function(m){
        var mm = ERA5_daily.filter(ee.Filter.calendarRange(m,m,'month'))
                           .select('total_precipitation')
                           .map(bigger10)
                           .map(countD10).sum()
                           .set('month',m).rename('D10');
        return mm;
    }));

    //Compute the Difference Per Month
    var R10 = bymonthR10.map(function(m){return m.multiply(ee.Number(7.5))});
    var D10 = bymonthD10.map(function(m){return m.multiply(ee.Number(150))});

    var lisR10  = R10.toList(R10.size());
    var lisD10  = D10.toList(D10.size());
    var seq = ee.List.sequence(0,lisR10.size().subtract(1));
    var diff = ee.ImageCollection(seq.map(function(i){return ee.Image(lisR10.get(i)).subtract(lisD10.get(i))}));
    
    //Multiply the Difference Per Month with the number of total erosive events which is equal to the D10 number
    var listbymonthD10 = bymonthD10.toList(bymonthD10.size());
    var listdiff= diff.toList(diff.size());
    var multipl = ee.ImageCollection(seq.map(function(i){return ee.Image(listdiff.get(i)).multiply(ee.Image(listbymonthD10.get(i)))}));
    
  return multipl;
}

//Annual R-factor for each of 30-years period, after checking that is positive
var R1989 = calc('1989-01-01','1990-01-01',box).map(positive).sum();
var R1990 = calc('1990-01-01','1991-01-01',box).map(positive).sum();
var R1991 = calc('1991-01-01','1992-01-01',box).map(positive).sum();
var R1992 = calc('1992-01-01','1993-01-01',box).map(positive).sum();
var R1993 = calc('1993-01-01','1994-01-01',box).map(positive).sum();
var R1994 = calc('1994-01-01','1995-01-01',box).map(positive).sum();
var R1995 = calc('1995-01-01','1996-01-01',box).map(positive).sum();
var R1996 = calc('1996-01-01','1997-01-01',box).map(positive).sum();
var R1997 = calc('1997-01-01','1998-01-01',box).map(positive).sum();
var R1998 = calc('1998-01-01','1999-01-01',box).map(positive).sum();
var R1999 = calc('1999-01-01','2000-01-01',box).map(positive).sum();
var R2000 = calc('2000-01-01','2001-01-01',box).map(positive).sum();
var R2001 = calc('2001-01-01','2002-01-01',box).map(positive).sum();
var R2002 = calc('2002-01-01','2003-01-01',box).map(positive).sum();
var R2003 = calc('2003-01-01','2004-01-01',box).map(positive).sum();
var R2004 = calc('2004-01-01','2005-01-01',box).map(positive).sum();
var R2005 = calc('2005-01-01','2006-01-01',box).map(positive).sum();
var R2006 = calc('2006-01-01','2007-01-01',box).map(positive).sum();
var R2007 = calc('2007-01-01','2008-01-01',box).map(positive).sum();
var R2008 = calc('2008-01-01','2009-01-01',box).map(positive).sum();
var R2009 = calc('2009-01-01','2010-01-01',box).map(positive).sum();
var R2010 = calc('2010-01-01','2011-01-01',box).map(positive).sum();
var R2011 = calc('2011-01-01','2012-01-01',box).map(positive).sum();
var R2012 = calc('2012-01-01','2013-01-01',box).map(positive).sum();
var R2013 = calc('2013-01-01','2014-01-01',box).map(positive).sum();
var R2014 = calc('2014-01-01','2015-01-01',box).map(positive).sum();
var R2015 = calc('2015-01-01','2016-01-01',box).map(positive).sum();
var R2016 = calc('2016-01-01','2017-01-01',box).map(positive).sum();
var R2017 = calc('2017-01-01','2018-01-01',box).map(positive).sum();
var R2018 = calc('2018-01-01','2019-01-01',box).map(positive).sum();
var R2019 = calc('2019-01-01','2020-01-01',box).map(positive).sum();

var all = R1989.add(R1990).add(R1991).add(R1992).add(R1993).add(R1994).add(R1995).add(R1996).add(R1997).add(R1998).add(R1999)
          .add(R2000).add(R2001).add(R2002).add(R2003).add(R2004).add(R2005).add(R2006).add(R2007).add(R2008)
          .add(R2009).add(R2010).add(R2011).add(R2012).add(R2013).add(R2014).add(R2015).add(R2016).add(R2017)
          .add(R2018).add(R2019);

var meanall = all.divide(30).rename('mean_R');

//LS-factor as an "Asset"
var lsfactor = ee.Image("users/alexandridisvasileios/ls_fac_cre_new").clip(box.geometry());
//K-factor as an "Asset"
var kfactor =  ee.Image("users/alexandridisvasileios/kfactor_crete").clip(box.geometry());
//P-factor as an "Asset"
var pfactor =  ee.Image("users/alexandridisvasileios/P_factor_crete_wgs84_fill").clip(box.geometry());

//Final A-factor is computed based on the 5 factors
var Afactor= meanall.multiply(cfactor_ndvi).multiply(lsfactor).multiply(kfactor).multiply(pfactor);

//----- Display Things&Layers----//
//Visualization
var viscfactor = {min:-1,max:1,palette:['green','white','black']}
var visrfactor = {min:0,max:2000,palette:['yellow','orange','purple']}
var visafactor = {min:0,max:500,palette:['yellow','orange','purple']}

Map.centerObject(box,8.3);
Map.addLayer(box.draw({color: '999999', strokeWidth: 2}),{},'Study Area');
Map.addLayer(Afactor,visafactor,'Afactor');

//Compute & Visualize some statistics values: min, max, mean
var min_value = Afactor.reduceRegion({
  reducer: ee.Reducer.min(),
  geometry: box.geometry(),
  scale: 10000,
  maxPixels: 1e16,
  // tileScale: 16
});
print(min_value);

var max_value = Afactor.reduceRegion({
  reducer: ee.Reducer.max(),
  geometry: box.geometry(),
  scale: 10000,
  maxPixels: 1e16,
  // tileScale: 16
});
print(max_value);

var mean_value = Afactor.reduceRegion({
  reducer: ee.Reducer.mean(),
  geometry: box.geometry(),
  scale: 10000,
  maxPixels: 1e16,
  // tileScale: 16
});
print(mean_value);

//----------------EXPORT----------------------------------
// Export A_parameter;
Export.image.toDrive({
  image: Afactor,
  description: 'Afactor',
  scale:100,
  maxPixels: 3784216672400,
});
