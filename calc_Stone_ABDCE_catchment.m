function [A,B,C,D,E,M] = calc_Stone_ABDCE_catchment(DEM,utmzone,Lat,a,b,c,d,e,m)

% get extent of raster in UTM
x_raster = DEM.georef.SpatialRef.XWorldLimits;
y_raster = DEM.georef.SpatialRef.YWorldLimits;

% get extent of raster in Lat Lon
[raster_lat,raster_lon]=utm2ll(x_raster,y_raster,utmzone,'wgs84');

% Create x and y array for interpolation of production rate within drainage
% basin in UTM coords
x_utm = DEM.refmat(3,1): DEM.cellsize: DEM.refmat(3,1)+ DEM.cellsize*(DEM.size(2)-1);
y_utm = DEM.refmat(3,2) - DEM.cellsize*(DEM.size(1)-1): DEM.cellsize : DEM.refmat(3,2);
y_utm = y_utm';

% Interpolation of scaling fator - to find values at [x_utm(:),y_utm(:)]
A = GRIDobj(x_utm,y_utm,repmat(interp1(y_raster,interp1(Lat,a,raster_lat),y_utm),1,length(x_utm)));
B = GRIDobj(x_utm,y_utm,repmat(interp1(y_raster,interp1(Lat,b,raster_lat),y_utm),1,length(x_utm)));
C = GRIDobj(x_utm,y_utm,repmat(interp1(y_raster,interp1(Lat,c,raster_lat),y_utm),1,length(x_utm)));
D = GRIDobj(x_utm,y_utm,repmat(interp1(y_raster,interp1(Lat,d,raster_lat),y_utm),1,length(x_utm)));
E = GRIDobj(x_utm,y_utm,repmat(interp1(y_raster,interp1(Lat,e,raster_lat),y_utm),1,length(x_utm)));
M = GRIDobj(x_utm,y_utm,repmat(interp1(y_raster,interp1(Lat,m,raster_lat),y_utm),1,length(x_utm)));
end

