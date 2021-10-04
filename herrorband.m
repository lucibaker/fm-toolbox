function [hl,hr] = herrorband(x,y,l,r,bandstyle,bandcolor)
% adds left and right horizontal error bands to a plot
% x and y: original x and y data
% l: left error bar length
% r: right error bar length
% bandstyle: line style of error bounds
hold on;
hr = plot(x + r, y, bandstyle, 'color', bandcolor);
hl = plot(x - l, y, bandstyle, 'color', bandcolor);