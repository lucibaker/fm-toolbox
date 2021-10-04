function [vl,vr] = errorband(x,y,t,b,bandstyle,bandcolor)
% adds left and right horizontal error bands to a plot
% x and y: original x and y data
% t: top error bar length
% b: bottom error bar length
% bandstyle: line style of error bounds
hold on;
vr = plot(x, y + t, bandstyle, 'color', bandcolor);
vl = plot(x, y - b, bandstyle, 'color', bandcolor);