# Modified_RUSLE_model

A modified RUSLE model has been implemented on Google Earth Engine (GEE) platform.

The model is entirely written in Javascript API where the required datasets are retrieved either from the GEE catalogue or uploaded as assets. In particular, satellite images Sentinel-2A and ERA-5 percipitation dataset are available in the data catalogue, and therefore the C and R factors were calculated, respectively. As for the K, LS and P factors were estimated on a GIS system and uploaded as assets in the platform. The K factor was obtained from the soil map of Greece, the LS from ALOS DEM and the P factor was directly derived in raster format from the ESDAC dataset.

This model is available after request, transferable to any region and applicable with other datasets as well (e.g. satellite images, precipitation data).
