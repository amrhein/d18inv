function [D tD Di core_order_new] = mkDni(path1,lims)

core_order_new = {
    'SA1967'
    'NA3223'
    'NA1299'
    'SA3770'

    'EI2100'
    'NA3146'    
    'NI1580'
    'EP3210'}';
files = dir([path1 '*.mat']);

lf = length(files);

core_order = {};
D = [];
tD = [];
Di = [];

for i = 1:lf
    
    % load the data and look at its properties
    eval(['load ' path1 files(i).name ';']);
    name = strrep(files(i).name,'.mat','');
    eval(['lat=' name '.lat' ';'])
    eval(['lon=' name '.lon' ';'])
    eval(['depth=' name '.depth' ';'])
    eval(['data=' name '.data' ';'])
    
    core_order = [core_order,name];
    
    % determine which age model will be used. If there are data in the
    % sixth column of the 'data' field, then I have constructed an age
    % model to compensate for otherwise missing data and that age model
    % should be used. Otherwise use the author's age models.
    if isempty(data)
        continue
    end
    if ~all(isnan(data(:,6)))
        t = data(~isnan(data(:,5)),6);
    else
        t = data(~isnan(data(:,5)),2);
    end

    % extract the d18O data
    d18 = data(~isnan(data(:,5)),5);
    
    % eliminate points for which t or d18O data are missing
    killem = (isnan(d18)|isnan(t));
    d18(killem) = [];
    t(killem) = [];

    % eliminate multiple measurements at the same time (the result of
    % samples at the same depth) by using the average
    
    allreps = [];
    for i2 = 1:length(t)
        repinds = find(t==t(i2));
        d18(repinds) = mean(d18(repinds));
        if numel(repinds)>1, allreps = [allreps;repinds(2:end)]; end
    end
    d18(unique(allreps)) = [];
    t(unique(allreps)) = [];
    
    % convert t to actual bp age
    t = t*1000;
    
    % flip the record so that time proceeds in the correct direction (in
    % units of ybp, so always 
    t = -flipud(t(:));
    d18 = flipud(d18(:));
    
    % omit times outside of those specified in lims
    d18(t<lims(1)|t>lims(2)) = [];
    t(t<lims(1)|t>lims(2)) = [];
    
    
    % Concatenate this record onto the D vector
    D = [D;d18(:)];
    tD = [tD;t(:)];
    Di = [Di;ones(length(d18),1)*i];
  
end
















