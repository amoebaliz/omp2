function drswitch(a)
%DRSWITCH.M: utility for omp2, allows directory change during programm execution

i = 'y';
while i == 'y'
	disp('  ')
	disp([a ' your current work directory, which is'])
	disp(['  ' pwd])
	disp('  ')
	disp('This directory contains the following directories/folders and files:')
	dir
	disp('  ')
	i = 'n';
	incontrol = input('Do you want to change directory (y/n)? [n]  ','s');
	if length(incontrol) > 0 i = incontrol; end
	if i == 'y'
		chdir = '..';
		incontrol = input('Change to which directory: [..]  ','s');
		if length(incontrol) > 0 chdir = incontrol; end
		newdir = ['cd ' chdir];
		eval(newdir);
	end
end
