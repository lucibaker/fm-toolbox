% print all figs in a folder as pngs
d = dir;
for i = 1:length(d)
if contains(d(i).name,'.fig')
openfig(d(i).name);
pause(1/10);
print(d(i).name(1:end-4),'-dpng','-r600');
end
end