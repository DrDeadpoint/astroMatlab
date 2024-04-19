hFig = figure(); hold on;
hAx = gca;
plot3DBody('earth', 0.3, [0,0,0]);
plot3DBody('moon', 0.2, [1.5, 0, 0]);
plot3DBody('sun', 0.8, [4, 3.7, 1.2]);
set(hAx, 'color', 'k', 'xTickLabels', {}, 'yTickLabels', {});
set(hFig, 'color', 'k');
axis equal;
axis off;
view(-57,0);