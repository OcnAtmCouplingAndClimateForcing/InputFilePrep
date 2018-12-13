%Average_Winter_Summer.m
%Emily Norton
%June 8, 2018
%This program opens up a netcdf file and calculates point-by-point averages over the the months and years of your choosing.  Note: Year data corresponds to the final month in the range. Averages from each year are then stored in a horizontal vector (ordered by lon), and all vectors are horizontally organized so each year has its own row.  Ultimately this matrix is exported as a csv (in addition to being saved as a .mat).  Also, a land mask is used to remove land-points.

clear all
close all

%Change the following variables
VAR = 'skt';  %'slp' = Sea Level Pressure; 'skt' = Skin Temp; 'hgt' = Geopotential Height at 200mbar; 'UFLX' or 'VFLX' for U/V wind stress; 'MLD' for mixed layer depth
YEARMIN = 1948;  %start year: summer(May-Sep,YEARMIN, avg), winter(Nov, YEARMIN-1 -Mar, YEARMIN); Note: MLD runs from 1959-2011
YEARMAX = 2018;  %end year
MOMIN = 11 ;  %start month - numerical input required; 'winter' is Nov - Mar (i.e. MOMIN = 11; MOMAX = 3); 'summer' is May - Sep (5 - 9); other periods of interest: slp from Nov - Jan, skt from Feb - Apr
MOMAX = 3;   %end month - numerical input required

% Select lat/lon bounds
latmin = 20;
latmax = 65;
lonmin = 120; 
lonmax = 255;

%%%-----------Shouldn't need to change below this line--------------
switch VAR
case 'skt'
FILENM = sprintf('/home/disk/clark/emilyln/FATE_Hunsicker_Bond_OcnAtmCoupling/RawData/SkinTemp_NCEP_%ito%i.nc', YEARMIN, YEARMAX);;

case 'slp'
FILENM = sprintf('/home/disk/clark/emilyln/FATE_Hunsicker_Bond_OcnAtmCoupling/RawData/SLP_NCEP_%ito%i.nc', YEARMIN, YEARMAX);

case 'hgt'
FILENM = sprintf('/home/disk/clark/emilyln/FATE_Hunsicker_Bond_OcnAtmCoupling/RawData/GeoPotHt_NCEP_%ito%i_200mbar.nc', YEARMIN, YEARMAX);

case 'UFLX'
FILENM = '/home/disk/clark/emilyln/FATE_Hunsicker_Bond_OcnAtmCoupling/RawData/U_Stress_NCEP_Full_Renamed.nc';
%FILENM = '/home/disk/clark/emilyln/FATE_Hunsicker_Bond_OcnAtmCoupling/RawData/U_Stress_NCEP_Renamed.nc';

case 'VFLX'
FILENM = '/home/disk/clark/emilyln/FATE_Hunsicker_Bond_OcnAtmCoupling/RawData/V_Stress_NCEP_Full_Renamed.nc';
%FILENM = '/home/disk/clark/emilyln/FATE_Hunsicker_Bond_OcnAtmCoupling/RawData/V_Stress_NCEP_Renamed.nc';

case 'MLD'
FILENM = '/home/disk/clark/emilyln/FATE_Hunsicker_Bond_OcnAtmCoupling/RawData/MLD_ORAS3_Renamed.nc';


end   %switch VAR filename

%%Open up netcdf file
ncid = netcdf.open(FILENM);
latid = netcdf.inqVarID(ncid,'lat');
lonid = netcdf.inqVarID(ncid,'lon');
timeid = netcdf.inqVarID(ncid,'time');
varid = netcdf.inqVarID(ncid,VAR);


%Extract Variables values from file
lat = netcdf.getVar(ncid,latid);
lon = netcdf.getVar(ncid,lonid);
time = netcdf.getVar(ncid,timeid);
var = netcdf.getVar(ncid,varid);

netcdf.close(ncid)

%--------------------------------
switch VAR
case 'skt'
%Apply detrend
detrendval = csvread('detrend_correction_skt.csv');
vardet = var;
for d = 1:length(detrendval)
	vardet(:,:,d) = var(:,:,d)-detrendval(d);
end	
var = vardet;
end


% Limit boundaries, and re-shape grids
varOrig = var;

% Re-shape grid, to match latG/lonG
varRe = permute(var,[2 1 3]);
var = varRe;

if lat(1)>lat(end)
[lonG, latG] = meshgrid(lon,lat(end:-1:1));
var = var(end:-1:1,:,:);
disp('lat is regular but not plaid. it has been re-ordered')
else 
[lonG, latG] = meshgrid(lon,lat);
disp('lat is regular and plaid')
end


latinds = find(latG <= latmax & latG >= latmin);
loninds = find(lonG <= lonmax & lonG >= lonmin);

inboundinds = intersect(latinds,loninds);
allinds = [1:length(latG(:))];
oobinds = setdiff(allinds,inboundinds);


%%Convert time to a meaningful date (currently saved as hours since 1 Jan 1800)

switch VAR 
case {'skt','slp'}
dayssinceref = time./24;
mattime = datenum('1 Jan 1800') + datenum(dayssinceref);  %for ref time = 1 Jan 1800
dataminyr = datevec(mattime(1));  %grab just the minimum year from the dataset
datamaxyr = datevec(mattime(end));  %grab just the max year from the dataset

case {'UFLX','VFLX','MLD'}
mattime = datenum('1/1/1') + datenum(time);   %this is what is indicated in the ncdump. even though it doesn't quite replicate the time signal, it is close enough
dataminyr = datevec(mattime(1));  
datamaxyr = datevec(mattime(end));

end  %switch VAR for date/time calibration

%For each point in the var matrix, average across summer (May Year - Sep,Year) and winter (Nov Year - Mar Year+1)
%first, pre-allocate final mats to save data in
numpointsAll = length(lat)*length(lon);
yrrange = YEARMIN:YEARMAX;
datMat = nan(length(yrrange),numpointsAll);

months = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
daymaxes = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];   %number of days in each month (excluding leap years)
mospan = MOMAX - MOMIN + 1;  %how many months we're averaging over
if mospan <0
mospan = mospan + 12;
end

for t = 1:length(yrrange)   %calculate summer and winter averages for each year
curryr = yrrange(t);

if MOMIN > MOMAX 
datsstr = sprintf('%i/%i/1',curryr-1,MOMIN);
else 
datsstr = sprintf('%i/%i/1',curryr,MOMIN);
end 
datestr = sprintf('%i/%i/%i',curryr,MOMAX,daymaxes(MOMAX));

datsdate = datenum(datsstr);
datedate = datenum(datestr);

datdateinds = find(datsdate<=mattime & datedate>=mattime);

datmat = var(:,:,datdateinds);

datavg = nanmean(datmat,3);

if size(datmat,3)~=mospan   %but if they didn't have the full complement of monthly data available, make them nans 
datavg = nan(size(datavg));
end

datMat(t,1:numpointsAll) = datavg(:);   %Note: The first or last row of datMat may be nans, depending on the range of months and the range of slp/skt data used 

%clear sumdateinds windateinds summat winmat sumavg winavg
end 


%--------------------------------
%%Identify land mask points and remove them from the matrices 

%%Open up land mask netcdf file
lmfile = '/home/disk/clark/emilyln/FATE_Hunsicker_Bond_OcnAtmCoupling/RawData/LandSeaMask_NCEP_Jan_15-70N_120-255E.nc'; 
ncidlm = netcdf.open(lmfile);
latidlm = netcdf.inqVarID(ncidlm,'lat');
lonidlm = netcdf.inqVarID(ncidlm,'lon');
timeidlm = netcdf.inqVarID(ncidlm,'time');
varidlm = netcdf.inqVarID(ncidlm,'lsmask');


%Extract Variables values from land mask file
latlm = netcdf.getVar(ncidlm,latidlm);
lonlm = netcdf.getVar(ncidlm,lonidlm);
timelm = netcdf.getVar(ncidlm,timeidlm);
varlm = netcdf.getVar(ncidlm,varidlm);  %In land mask, sea = 0, land = -1

netcdf.close(ncidlm)

%Since the land mask grid is different from the SLP and SST grid, we'll need to do some interpolation to mask out the land.

if size(varlm) == size(squeeze(var(:,:,1)))
landma = varlm;
else
[lonlmG,latlmG] = meshgrid(lonlm,latlm);
%[lonG,latG] = meshgrid(lon,lat);
landma = interp2(lonlmG,latlmG,varlm',lonG,latG);
end
seainds = find(landma>=-0.5 & landma<=0);  %We might want to change this threshold for what qualifies as sea... most stringent would be landmask ==0
%landinds = find(landma == -1);
landinds = find(landma~=0);


% Also identify any flag values
switch VAR
case 'MLD'
nonflaginds = find(datMat(60,:)<1000000);     %will need to vary the row, needs to have data present; NOTE: this assumes that the flgged inds are the same in each row...which should be ~true for a land mask, but won't be true for a physical/physical data set, e.g.

otherwise 
nonflaginds = allinds;
end  %switch VAR for nanning out flag values

%Now grab only those seainds, inboundinds, and (if applicable) non-flag inds from datMat

goodinds = intersect(seainds,inboundinds);
goodinds = intersect(goodinds, nonflaginds);

datSeaMatNan = datMat;   %use this as a test: keep dimensions the same as the grid, just nan out out of bounds inds and land inds; just note that the kohonen package in R doesn't accept entire columns of nans

for t = 1:length(datMat(:,1))
%datSeaMat(t,1:length(seainds)) = datMat(t,seainds);
datSeaMat(t,1:length(goodinds)) = datMat(t,goodinds);
datSeaMatNan(t,landinds) = nan;   			%basically doing the same thing, but keeping the full grid size
datSeaMatNan(t,oobinds) = nan;
end


% Check that the land mask and the variable field looks good
testSeaMatNan = reshape(datSeaMatNan(60,:),size(var(:,:,1)));

figure
contourf(lonG,latG,testSeaMatNan)
title('testSeaMatNan')
colorbar



%--------------------------------------
%%Save output as .mat and csv
datfilenmsav = sprintf('MonthlyAvg%ito%i_forYears%ito%i_Bounds%ito%iN_%ito%iE_%s_detrend_fewerNans_incLat.csv',MOMIN,MOMAX,YEARMIN,YEARMAX,latmin,latmax,lonmin,lonmax,VAR);
datafilenmsav = sprintf('MonthlyAvg%ito%i_forYears%ito%i_Bounds%ito%iN_%ito%iE_%s',MOMIN,MOMAX,YEARMIN,YEARMAX,latmin,latmax,lonmin,lonmax,VAR);

csvwrite(datfilenmsav,datSeaMat);
%csvwrite(datfilenmsav,datSeaMatNan);

%save(datafilenmsav)


% Make other lat/lon/seaind files for plotting later
lonBind = find(lon <= lonmax & lon >= lonmin);
lonB = lon(lonBind);
latBind = find(lat <= latmax & lat >= latmin);
latB = lat(latBind);

latcsv = sprintf('%s_lat_vec_Bounds%ito%iN_%ito%iE.csv', VAR, latmin, latmax, lonmin, lonmax);
%csvwrite(latcsv,latB);

loncsv = sprintf('%s_lon_vec_Bounds%ito%iN_%ito%iE.csv', VAR, latmin, latmax, lonmin, lonmax);
%csvwrite(loncsv,lonB);

seacsv = sprintf('%s_seainds_Bounds%ito%iN_%ito%iE.csv', VAR, latmin, latmax, lonmin, lonmax);
%csvwrite(seacsv,goodinds);


% if we're keeping the whole grid, and just nanning out the land mask and out of bounds inds
latcsv = sprintf('%s_lat_vec_Bounds_OrigGrid_20to65N_120to255E.csv', VAR);
%csvwrite(latcsv,lat);

loncsv = sprintf('%s_lon_vec_Bounds_OrigGrid_20to65N_120to255E.csv', VAR);
%csvwrite(loncsv,lon);

seacsv = sprintf('%s_seainds_Bounds_OrigGrid_%ito%iN_%ito%iE_incLat.csv', VAR, latmin, latmax, lonmin, lonmax);
csvwrite(seacsv,goodinds);

