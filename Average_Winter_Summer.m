%Average_Winter_Summer.m
%Emily Norton
%June 8, 2018
%This program opens up a netcdf file and calculates point-by-point averages over the summer (May-Sep) and winter (Nov-Mar) for a given range of years.  Note: Year data corresponds to the final month in the range. Averages from each year are then stored in a horizontal vector (ordered by lon), and all vectors are horizontally organized so each year has its own row.  Ultimately this matrix is exported as a csv (in addition to being saved as a .mat).  Also, a land mask is used to remove land-points (TBD when/where this happens).


%Change the following variables
VAR = 'slp';  %'slp' = Sea Level Pressure; 'skt' = Skin Temp
YEARMIN = 1989;  %start year: summer(May-Sep,YEARMIN, avg), winter(Nov, YEARMIN-1 -Mar, YEARMIN)
YEARMAX = 2018;  %end year

%FILENM = '/home/disk/clark/emilyln/FATE_Hunsicker_Bond_OcnAtmCoupling/SkinTemp_NCEP_1948to2018.nc';
FILENM = '/home/disk/clark/emilyln/FATE_Hunsicker_Bond_OcnAtmCoupling/SLP_NCEP_1948to2018.nc';

%%%-----------Shouldn't need to change below this line--------------

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
sumMat = nan(length(yrrange),numpointsAll);
winMat = nan(length(yrrange),numpointsAll);


for t = 1:length(yrrange)   %calculate summer and winter averages for each year
curryr = yrrange(t);

sumsstr = sprintf('1 May %i',curryr);
sumestr = sprintf('30 Sep %i',curryr);
winsstr = sprintf('1 Nov %i',curryr-1);
winestr = sprintf('31 Mar %i',curryr);

sumsdate = datenum(sumsstr);
sumedate = datenum(sumestr);
winsdate = datenum(winsstr);
winedate = datenum(winestr);

sumdateinds = find(sumsdate<=mattime & sumedate>=mattime);
windateinds = find(winsdate<=mattime & winedate>=mattime);

summat = var(:,:,sumdateinds);
winmat = var(:,:,windateinds);

sumavg = nanmean(summat,3);
winavg = nanmean(winmat,3);

if size(summat,3)~=5   %but if they didn't have the full complement of monthly data available, make them nans 
sumavg = nan(size(sumavg));
end
if size(winmat,3) ~=5
winavg = nan(size(winavg));
end

sumMat(t,1:numpointsAll) = sumavg(:);   %Note: The last row of sumMat will be nans, because summer not calculated for last year in the range
winMat(t,1:numpointsAll) = winavg(:);  %Note: The first row of winMat will be nans so as to keep the year consistent with sumMat (with year pertaining to the last month in the seasonal range)

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

%Now grab only those seainds from sumMat and winMat

for t = 1:length(sumMat(:,1))
sumSeaMat(t,1:length(seainds)) = sumMat(t,seainds);
winSeaMat(t,1:length(seainds)) = winMat(t,seainds);
end


%--------------------------------------
%%Save output as .mat and csv
sumfilenmsav = sprintf('SummerAvg_%ito%i_%s_TEST2.csv',YEARMIN,YEARMAX,VAR);
winfilenmsav = sprintf('WinterAvg_%ito%i_%s_TEST2.csv',YEARMIN,YEARMAX,VAR); 
datafilenmsav = sprintf('WinterSummerAvg_%ito%i_%s_TEST2',YEARMIN,YEARMAX,VAR);

csvwrite(winfilenmsav,winSeaMat);
csvwrite(sumfilenmsav,sumSeaMat);

save(datafilenmsav)



