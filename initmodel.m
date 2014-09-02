function model = initmodel(model);

	% Initialize the musculoskeletal model

	% initialize the model parameters
	warning('OFF','MATLAB:xlsread:Mode');
	model.gait2d('Initialize', xlsread(model.parameterfile,'','','basic') );
	warning('ON','MATLAB:xlsread:Mode');
	
	% if requested, amputate the model
	if strcmp(model.type, 'bka')
		model.gait2d('Set','Right BKA', model.anklestiffness);		
	elseif strcmp(model.type, 'aka')
		model.gait2d('Set','Right AKA', model.anklestiffness);	
	end
	
	% ask what the mass is, and store it
	model.mass = model.gait2d('Get','Total Mass');
	% fprintf('Total mass of model: %8.3f kg\n',model.mass);
	
	% ask for the Fmax values of the muscles
	model.Fmax = model.gait2d('Get','Fmax');

end
