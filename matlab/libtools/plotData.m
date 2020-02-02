function  plotData(ndat,data,caption,leg,sfxData,opt)
% Plot multiple data sets given in the argument as 2-column array X,Y. 
% A 3rd column with errors is optional.
% The data can contain multiple sets, partitioning into data sets is defined
% by the ndat array.
%
% opt = Options structure
% Option.groups = array of integers. Number of elements >= number of data sets.
% Each value defines corresponding color, point style and line style, 
% which are predefined in the colors, styles and ptstyles arrays.
%
% Option.types: similarly to groups, types is an array of integers which
% defines curve style, 0=line, 1=points, 2=points and line
%
% leg = set of legend items
%
% sfxData = string added to each legend item
%
% opt = Options structure
%
% version: 1.3
% date: 2/8/2016
%

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1,'YGrid',opt.YGrid,'XGrid',opt.XGrid,...
    'YMinorTick',opt.YMinorTick,'XMinorTick',opt.XMinorTick);
% Set limits if required
if (opt.xauto ~= true); xlim(axes1,opt.xlimits);end;
if (opt.yauto ~= true); ylim(axes1,opt.ylimits);end;
box(axes1,'on');
hold(axes1,'all');

ny=numel(ndat);
ncol=size(data,2);
X=data(:,1);
Y=data(:,2);
if (ncol>2)
    YE=data(:,3);
else
    YE=0;
end;
% plot data
i1=1;
for i=1:ny
  i2=i1-1+ndat(i);  
  igr=opt.groups(i);
  col=opt.colors(igr,:);
  cap=[char(leg(1,i)) ' ' sfxData];
  ityp=opt.types(i); 
  if (ncol>2)
      plt1=errorbar(X(i1:i2),Y(i1:i2),YE(i1:i2),'Parent',axes1,'Color',col,'DisplayName',cap);
  else
      plt1=plot(X(i1:i2),Y(i1:i2),'Parent',axes1,'Color',col,'DisplayName',cap);
  end;
  if (ityp==0)        
     sty=char(opt.styles(igr));
     set(plt1,'LineStyle',sty,'LineWidth',opt.LineWidth);
  elseif (ityp==1)
     sty=char(opt.ptstyles(igr));
     set(plt1,'Marker',sty,'LineStyle','none');
  else
     sty=char(opt.ptstyles(igr));
     set(plt1,'Marker',sty,'LineStyle','-','LineWidth',1);  
  end;
  i1=i1+ndat(i);
end

% Create labels
xlabel(opt.xlabel);
ylabel(opt.ylabel);

% Create title
title({caption});

% Create legend
if ( size(char(leg),2)>0)
legend1 = legend(axes1,'show');
set(legend1,'Location',opt.LegendPosition,'FontSize',10);
end


