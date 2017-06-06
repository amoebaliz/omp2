function ODV2OMP2(infile,outfile)
%
% Usage ODV2OMP2(infile, outfile)
% 
% The routine ODV2OMP1 reads a text file produced by ODV in the ODV
% spreadsheet format, and produces an output mat file with the name
% 'outfile', in the form required for input to the OMP routines.
% This program keeps only the values that contain all the biochemical
% parameters (dissolved oxygen, silicates, nitrates and phosphates)
%
% Vassilis Zervakis and Stallo Leontiou, Mytilene, University of the
% Aegean, 2010
%
fid = fopen(infile);
strline = fgetl(fid);
hydrodata = []
while strline ~= -1,
  if strline(1:2) == '//'
      strline = fgetl(fid);
  else
      if strline(1:6) == 'Cruise',
          strline = fgetl(fid);
      end
      if strline(1) == '	',
          % find where in the line there are tabs
          jtabs = [];
          jtabs = find(strline ~= '	');
          % find the parts where non-tabs are continuous (so they contain
          % data)
          ilast = max(jtabs);
          istart = min(jtabs);
          ivalues = 0;
          i = 1;		
          ibegin = i;
          numtimes = 0;
          while(jtabs(i) < ilast),
             ivalues = ivalues + 1;
             if jtabs(i)+1 == jtabs(i+1),
                 i = i+1;
             else
                 numtimes = numtimes+1;
                 iend = i;
                 value = strline(jtabs(ibegin):jtabs(iend));
                 i = i + 1;
                 ibegin = i;
                 if numtimes == 1,
                     depth = str2num(value);
                 end
                 if numtimes == 3,
                     temp = str2num(value);
                 end
                 if numtimes == 5,
                     sal = str2num(value);
                 end
                 if numtimes == 7,
                     oxy = str2num(value);
                 end
                 if numtimes == 9,
                     si = str2num(value);
                 end
                 if numtimes == 11,
                     nit = str2num(value);
                 end
                 if numtimes == 13,
                     ph = str2num(value);
                     if ((isfinite(temp)) & (temp ~= 1.)),
                         if (isfinite(oxy))&(isfinite(si))&(isfinite(nit))&(isfinite(ph)),
                             dat = [lat long depth temp sal oxy si nit ph];
                             if (length(lat) > 1),
                                 lat
                             end
                             hydrodata = [hydrodata ; dat];
                         end
                     end
                     numtimes = 0;
                 end
             end
          end
      else
          % find where in the line there are tabs
          jtabs = [];
          jtabs = find(strline ~= '	');
          % find the parts where non-tabs are continuous (so they contain
          % data)
          ilast = max(jtabs);
          istart = min(jtabs);
          ivalues = 0;
          i = 1;		
          ibegin = i;
          numtimes = 0;
          jtabs;
          while(jtabs(i) < ilast),
             ivalues = ivalues + 1;
             if jtabs(i)+1 == jtabs(i+1),
                 i = i+1;
             else
                 numtimes = numtimes+1;
                 iend = i;
                 value = strline(jtabs(ibegin):jtabs(iend));
                 i = i + 1;
                 ibegin = i;
                 if numtimes == 5,
                     long = str2num(value);
                 end
                 if numtimes == 6,
                     lat = str2num(value);
                     A = [lat long]
                 end
                 end
             end
          end
      end
      strline = fgetl(fid);
  end
end
% save everything in the output file
lat = hydrodata(:,1);
long = hydrodata(:,2);
depth = hydrodata(:,3);
temp = hydrodata(:,4);
sal = hydrodata(:,5);
oxy = hydrodata(:,6);
si = hydrodata(:,7);
nit = hydrodata(:,8);
ph = hydrodata(:,9);
% outfile = strcat(outcode,'.mat')
eval (['save ' outfile ' lat long depth temp sal oxy si nit ph']);
return