# River networks from OpenStreetMap data

This code uses OpenStreetMaps geometries to find and correct places within deltas where the RGIS STN  river networks are incorrect.

The STN river network is derived from a digital elevation model and assumes that rivers only converge as they flow downhill toward a receiving basin. This flow-direction algorithm works very well in regions with subtaintial relief, but fails in the flat floodplains of river deltas. Within deltas you typically find that rivers bifurcate in difficult-to-predict places.

This code generates a list of bifurcation points where the STN network should be modifed to divert some fraction of discharge to a neighboring node. Discharge fractionation is a function of river width, estimated here based on area and perimeter of river geomoetries, as well as extracted from the Global River Width from Landsat dataset.

Some of the code may need to be edited to run elsewhere as some paths to required datasets are hardcoded.

Use `runall.sh`, or run a specific spatial resolution with `STNres=06min scons`
