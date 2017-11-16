function currday = convertDate(dates,startdate,dateformat)
%This helper function converts dates to numbers, and then re-scales them to
%whole numbers of days since startdate

if nargin<3
    dateformat = 'yyyy-mm-dd  HH:MM:SS';
end



 datenumbers = datenum(dates, dateformat);
 whole_dates = floor(datenumbers);
 
 if nargin<2
     startdate = min(whole_dates);
 end
 
 currday = whole_dates - startdate+1;