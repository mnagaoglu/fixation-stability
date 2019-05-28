function [isoa, bcea, PRL, PRL2, density, xGrid, yGrid, fh] = ...
    ComputeFixationStability(xDeg, yDeg, cumProb, isShowPlot)
%[isoa, bcea, density, stats, fh] = ...
%    ComputeFixationStability(xDeg, yDeg, cumProb, isShowPlot)
%
%
% MNA 5/26/2018 wrote it. mnagaoglu@gmail.com
% MNA 5/28/2019 better input handling.
%
if nargin<3 || isempty(cumProb)
    cumProb = 0.68;
end

if nargin<4 || isempty(isShowPlot)
    isShowPlot = 1;
end


xy = [xDeg yDeg];


% set a cutoff for the trimming function. This is needed when you know you
% have outlier data points. For instance, in video-based eye trackers,
% before and after blinks, there will be some spuriously larger eye
% position data. By using a trim threshold, these outliers can be removed.
% By defaul it's set to 0. If you want to remove the upper 1%, set it to
% 0.01.
kdetrim = 0.0; 

% bandwidth size
bw = std(xy)/length(xy)^(1/6); 

try
    % Perform 2d kernel density estimation on the gaze data
    [~,density,xGrid,yGrid] = kde2d_trimmed(kdetrim, xy, 256,[], [], bw);

catch errkde
    fprintf('Not sufficient data for kernel density estimation!!!!\n');
    rethrow(errkde)
end

if any(isnan(density))
    isoa = [];
    bcea = [];
    density = [];
    xGrid = [];
    yGrid = [];
    PRL = [NaN NaN];
    PRL2 = PRL;
    fh = [];
    fprintf('Not sufficient data for kernel density estimation!!!!\n');
    return;
end


% compute BCEA and PRL
PRL = bimean(xGrid,yGrid,density);

% use the following way instead of the old-school method for BCEA. This is
% because, if 'kdetrim' above is set to a value >0, some of the data will
% be trimmed. We would want to use the trimmed version for both BCEA and
% ISOA to have directly comparable results. But in any case, this should
% not be much different than the old-school method. For reference, you can
% see the old-school method below (commented out).
pv = bivar(xGrid,yGrid,density);
bcea = pi*chi2inv(cumProb,2)*sqrt(prod(pv)); 

% old school BCEA method
% bcea = pi*chi2inv(cumProb,2)*sqrt(1-corr(xDeg,yDeg)^2)*nanstd(xDeg)*nanstd(yDeg); 

% compute PRL based on max density
maxdensity = max(density(:));
I = find(maxdensity == density);
PRL2 = [xGrid(I), yGrid(I)];


% show plot if necessary 
if isShowPlot
    fh = figure('visible','on');
else
    fh = figure('visible','off');
end
msize = 20;
mcolor = 'k';
sh = scatter(xDeg, yDeg,msize,mcolor,'filled','MarkerEdgeColor','none');
hold on;
sh.MarkerFaceAlpha = 0.2;
set(gca,'fontsize',16);
xlabel('Hor. position (deg)')
ylabel('Ver. position (deg)')

% plot contour map for the density using N levels
N = 25;
steps = linspace(min(density(:)),max(density(:)),N+2);
steps = steps(2:end-1);
[C, ch] = contour(xGrid,yGrid,reshape(density,size(xGrid,1),size(xGrid,2)),steps); hold on
ph(1) = plot(PRL(1),PRL(2),'+r','MarkerSize',15);
ph(2) = plot(PRL2(1),PRL2(2),'+g','MarkerSize',15);
caxis([0 .075]);
colormap(jet)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the numerical isoline method:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% basically, estimate the area of each iso-contour line and then
% interpolate to the isoline whose cumulative probabilty will be equal to
% a desired value (e.g., 0.68)

sumDensity = sum(density(:));
isoAreas = nan(N,1);
isoProb = nan(N,1);

% go over each iso-probability line and compute the area within it.
for i=1:length(steps)

    try
        if i==1
            st = 2; %#ok<NASGU>
            en = C(2,1)+1;
        else
            st = en+2;
            en = C(2,st-1)+st-1;
        end
        
        islands = find(C(1,:) == steps(i));
        in = false(length(xGrid(:)),1);
        tempArea = 0;
        for j=1:length(islands)   
            st = islands(j)+1;
            xv = C(1,st:C(2,st-1)+st-1);
            yv = C(2,st:C(2,st-1)+st-1);
            in = in | inpolygon(xGrid(:),yGrid(:),xv',yv');
            tempArea = tempArea + polyarea(xv,yv);
        end

        isoProb(i) = sum(density(in))/sumDensity;
        isoAreas(i) = tempArea;
        
    catch err
        err.message
        err.stack.line
        err.stack.name
        fprintf('Error occured during isoline method\n');
        rethrow(err)
    end
end


% interpolate to desired isoline 
try
    desiredStep = interp1(isoProb, steps, cumProb, 'pchip');
    figure(fh);
    delete(ch);
    C = contour(xGrid,yGrid,reshape(density,size(xGrid,1),size(xGrid,2)),...
        [desiredStep 1],'LineWidth',3);
    colormap([0 0 1])


    islands = find(C(1,:) == desiredStep);
    isoa = 0;
    for i=1:length(islands)
        st = islands(i)+1;
        xv = C(1,st:C(2,st-1)+st-1);
        yv = C(2,st:C(2,st-1)+st-1);
        isoa = isoa + polyarea(xv,yv);
    end
    stats.isoProb = isoProb;
    stats.isoAreas = isoAreas;
    
catch erriso
    erriso.message
    fprintf('\n\nError during interpolation or isoline area computation...\n')
    rethrow(erriso)
end


% if plot is visible, show the quantitative results on the figure as well
if isShowPlot
    figure(fh);
    text(min(xDeg)-0.8,min(yDeg)+0.2,...
        sprintf('PRL_m_u %.2f %.2f \n PRL_m_a_x_(_p_d_f_) %.2f %.2f \n BCEA: %.2f, ISOA: %.2f',...
        PRL(1), PRL(2), PRL2(1), PRL2(2),bcea, isoa));
    title(sprintf('Selected cumulative probabilty: %.2f',cumProb));
    
    xlim([min(xDeg) max(xDeg)] + [-1 2]);
    ylim([min(yDeg) max(yDeg)] + [-1 2]);
    inset = axes(gcf,'Position',[0.7 0.7 0.2 0.2]);
    plot(isoProb,isoAreas,'o-k','LineWidth',2);
    xlabel('Cum. Prob.')
    ylabel('ISOA')
    set(gca,'fontsize',14)
    
    legend(ph,{'PRL_m_u','PRL_m_a_x_(_p_d_f_)'},'location','northwest')
else
    if ~isempty(fh)
        delete(fh);
    end
    fh = [];
end


