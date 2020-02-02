function myfigure(X1, Y1, caption, opt, leg)
%% Custom line plot
% Using predefined option structure on input. |opt.groups| defines
% references to data groups, so that multiple data sets (columns of |Y1|)
% can share the same color, line style, legend item etc. For example, if
% opt.groups(1)=2 and opt.groups(3)=2, then the data sets Y1(:,1) and
% Y1(:,3) will use the same color, opt.colors(2) and so on for legend
% items, line styles etc.
%
% Use Utils.plotOptions() to generate default options
%

% get number of datastes
sz=size(Y1);
ncol=sz(2);          % number of datasets
ngr=max(opt.groups); % number of data groups

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

% Create multiple lines using matrix input to plot
plot1 = plot(X1,Y1,'Parent',axes1);
hasLegend=zeros(ngr,1);
showLegend=0;
cap='';
if (nargin>4) 
    showLegend=1;
end
for i=1:ncol
    igr=opt.groups(i);
    ityp=opt.types(i);
    if (showLegend==1)
        cap=char(leg(1,igr));
    end;
    col=opt.colors(igr,:);
    if (ityp==0 || ityp==2)        
        sty=char(opt.styles(igr));
        set(plot1(i),'LineStyle',sty,'DisplayName',cap,'Color',col,'LineWidth',opt.LineWidth);
    end
    if (ityp==1 || ityp==2)        
        sty=char(opt.ptstyles(igr));
        set(plot1(i),'Marker',sty,'DisplayName',cap,'Color',col);
    end;
    if (hasLegend(igr) ==0)
        hasLegend(igr)=1*showLegend;
    else
        % Exclude line from legend
        set(get(get(plot1(i),'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); 
    end;
end;

% Create labels
xlabel(opt.xlabel);
ylabel(opt.ylabel);


% Create title
title({caption});

% Create legend
if (opt.showlegend*showLegend)
legend1 = legend(axes1,'show');
set(legend1,'Location',opt.LegendPosition,'FontSize',12);
end

% set(get(get(hLine,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend


