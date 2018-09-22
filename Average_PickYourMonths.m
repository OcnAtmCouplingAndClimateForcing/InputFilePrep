%Average_Winter_Summer.m
%Emily Norton
%June 8, 2018
%This program opens up a netcdf file and calculates point-by-point averages over the the months and years of your choosing.  Note: Year data corresponds to the final month in the range. Averages from each year are then stored in a horizontal vector (ordered by lon), and all vectors are horizontally organized so each year has its own row.  Ultimately this matrix is exported as a csv (in addition to being saved as a .mat).  Also, a land mask is used to remove land-points.

clear all
close all

%Change the following variables
VAR = 'skt';  %'slp' = Sea Level Pressure; 'skt' = Skin Temp
YEARMIN = 1948;  %start year: summer(May-Sep,YEARMIN, avg), winter(Nov, YEARMIN-1 -Mar, YEARMIN)
YEARMAX = 2018;  %end year
MOMIN = 2 ;  %start month - numerical input required; 'winter' is Nov - Mar; 'summer' is May - Sep; other periods of interest: slp from Nov - Jan, skt from Feb - Apr
MOMAX = 4;   %end month - numerical input required


%%%-----------Shouldn't need to change below this line--------------
switch VAR
case 'skt'
FILENM = '/home/disk/clark/emilyln/FATE_Hunsicker_Bond_OcnAtmCoupling/SkinTemp_NCEP_1948to2018.nc';

case 'slp'
FILENM = '/home/disk/clark/emilyln/FATE_Hunsicker_Bond_OcnAtmCoupling/SLP_NCEP_1948to2018.nc';
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
%%Convert time to a meaningful date (currently saved as hours since 1 Jan 1800)
dayssinceref = time./24;
mattime = datenum('1 Jan 1800') + datenum(dayssinceref);  %for ref time = 1 Jan 1800
dataminyr = datevec(mattime(1));  %grab just the minimum year from the dataset
datamaxyr = datevec(mattime(end));  %grab just the max year from the dataset

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
lmfile = 'LandSeaMask_NCEP_Jan.nc'; 
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
landmask = varlm;
else
[lonlmG,latlmG] = meshgrid(lonlm,latlm);
[lonG,latG] = meshgrid(lon,lat);
landma = interp2(lonlmG,latlmG,varlm',lonG,latG);
landmask = landma';
end
seainds = find(landmask>=-0.5 & landmask<=0);  %We might want to change this threshold for what qualifies as sea... most stringent would be landmask ==0

%Now grab only those seainds from datMat

for t = 1:length(datMat(:,1))
datSeaMat(t,1:length(seainds)) = datMat(t,seainds);
end


%--------------------------------------
%%Save output as .mat and csv
datfilenmsav = sprintf('MonthlyAvg%ito%i_forYears%ito%i_%s.csv',MOMIN,MOMAX,YEARMIN,YEARMAX,VAR);
datafilenmsav = sprintf('MonthlyAvg%ito%i_forYears%ito%i_%s',MOMIN,MOMAX,YEARMIN,YEARMAX,VAR);

csvwrite(datfilenmsav,datSeaMat);

save(datafilenmsav)



