function dp = defaultPlotting
dp = struct;
dp.lineColor = colour('k');
dp.lineWidth = 3;
dp.lineStyle = '-';
dp.HandleVisibility = 'on';
dp.low_thrust = struct(); %fill with details about plotting thrust lines
dp.markerStyle = '*';
dp.markerSize = 15;
dp.LTplotmult = 0.1;
dp.LTplotstep = 2;
dp.LTdir = true;
dp.LTcolor = colour('k');
end