colmap = hsv(1600);

clf
for ii = 0:1:1599
    DELTA = readmatrix("bins_1600_rN_16/impulse_functions/delta_" + ii + ".txt");
    x = (DELTA(:,1)-DELTA(1,1)) / (DELTA(end,1)-DELTA(1,1));
    y = DELTA(:,3);

    p = polyfit(x, y, 4);

    x1=linspace(x(1),x(end));
    y1=polyval(p,x1);
    
    plot(x1, y1, 'Color', colmap(ii+1,:))
    ylim([-0.1 0.1])
    hold on

    frame = getframe(gcf);
    img = frame2im(frame);
    [img,cmap] = rgb2ind(img,256);
    %if ii == 0
    %    imwrite(img,cmap,'kick_[-0.1,0.1].gif','gif','LoopCount',Inf,'DelayTime',0.01);
    %else
    %    imwrite(img,cmap,'kick_[-0.1,0.1].gif','gif','WriteMode','append','DelayTime',0.01);
    %end

end
hold off