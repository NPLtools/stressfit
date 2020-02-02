function  plotFit(ncY,X,Y,ncF,XF,YF,caption,leg,sfxData,sfxFit,opt)
%PLOTFIT Plot multiple lines, designed for comparing data and fit
% Data are arrays with possibly multiple columns. If number of columns in X
% equals the number of columns in Y, each curve will use its own X-scale.
% If X has only 1 column, it will be used for all Y columns.
%
% Corresponding data and fit pairs use the same color, but different
% styles. Styles are defined by strings. Colors are given in opt.colors.
% ncom = partitioning into datasets: all datasets are in single X,Y,...
% columns.
%
% version: 1.2
% date: 4/4 2014
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

%nx=size(X,2);
ny=numel(ncY);
%nxf=size(XF,2);
%nyf=size(YF,2);

% plot data
i1=1;
for i=1:ny
  i2=i1-1+ncY(i);  
  igr=opt.groups(i);
  col=opt.colors(igr,:);
  cap=[char(leg(1,igr)) ' ' sfxData];
  ityp=opt.types(i);  
  plt1=plot(X(i1:i2),Y(i1:i2),'Parent',axes1,'Color',col,'DisplayName',cap);
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
  i1=i1+ncY(i);
end

% plot fit
ny=numel(ncF);
i1=1;
for i=1:ny
  i2=i1-1+ncF(i); 
  igr=opt.groups(i);
  col=opt.colors(igr,:);
  if (size(sfxFit,2)>0)
    cap=[char(leg(1,igr)) '  '  sfxFit];
  else
    cap='';
  end;  
  plt1=plot(XF(i1:i2),YF(i1:i2),'Parent',axes1,'Color',col,'DisplayName',cap);
  %sty=char(opt.styles(igr));
  set(plt1,'LineStyle','--','LineWidth',opt.LineWidth);  
  if (size(sfxFit,2)<=0)
  % Exclude line from legend
     set(get(get(plt1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); 
  end
  i1=i1+ncY(i);
end

% Create labels
xlabel(opt.xlabel);
ylabel(opt.ylabel);

% Create title
title({caption});

% Create legend
legend1 = legend(axes1,'show');
set(legend1,'Location',opt.LegendPosition,'FontSize',10);


