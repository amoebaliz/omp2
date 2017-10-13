function ax = displaytext(ht)

global ftm ftz fts;

ax = axes('position',[0.0 0.17 1 0.25],'color',[0.8 0.9 0.9], ...
	'xcolor',[0 0 0],'ycolor',[0 0 0],'xtick',[],'ytick',[]);
cla;
for i = 1:8
	text(0.01,1.005-0.115*i,ht(i,:),'color','k','fontname',eval(ftm),'fontsize',ftz);
end
