function  plotFitErr( X,Y,XF,YF,EF,caption,leg,sfxData,sfxFit,opt)
%PLOTFITERR Plot multiple lines, designed for comparing data and fit
% Data are arrays with possibly multiple columns. If number of columns in X
% equals the number of columns in Y, each curve will use its one X-scale.
% If X has only 1 column, it will be used for all Y columns.
%
% Corresponding data and fit pairs use the same color, but different
% styles. Styles are defined by strings. Colors are given in opt.colors.

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

nx=size(X,2);
ny=size(Y,2);
nxf=size(XF,2);
%nyf=size(YF,2);

% plot data
for i=1:ny
  ix=min(nx,i);
  igr=opt.groups(i);
  col=opt.colors(igr,:);
  cap=[char(leg(1,igr)) ' ' sfxData];
  ityp=opt.types(i);
  plt1=plot(X(:,ix),Y(:,i),'Parent',axes1,'Color',col,'DisplayName',cap);
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
end

% plot fit
for i=1:ny
  ix=min(nxf,i);
  igr=opt.groups(i);
  col=opt.colors(igr,:);
  if (size(sfxFit,2)>0)
    cap=[char(leg(1,igr)) '  '  sfxFit];
  else
    cap='';
  end;  
  plt1=errorbar(XF(:,ix),YF(:,i),EF(:,i),'Parent',axes1,'Color',col,'DisplayName',cap);
  %sty=char(opt.styles(igr));
  set(plt1,'LineStyle','--','LineWidth',opt.LineWidth);  
  if (size(sfxFit,2)<=0)
  % Exclude line from legend
     set(get(get(plt1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); 
  end
end

% Create labels
xlabel(opt.xlabel);
ylabel(opt.ylabel);

% Create title
title({caption});

% Create legend
legend1 = legend(axes1,'show');
set(legend1,'Location',opt.LegendPosition,'FontSize',10);


