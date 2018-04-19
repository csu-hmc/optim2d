function model = initmodel(model);

	% Initialize the musculoskeletal model

	% initialize the model parameters
	warning('OFF','MATLAB:xlsread:Mode');
    params = xlsread(model.parameterfile,'','','basic');
    if strcmp(model.type, 'torque')
        params(27+(1:16),2) = 0; %remove muscles
    end
	model.gait2d('Initialize', params );
	warning('ON','MATLAB:xlsread:Mode');
	
	% if requested, amputate the model
	if strcmp(model.type, 'bka')
		model.gait2d('Set','Right BKA', model.anklestiffness, model.ankledamping);
    elseif strcmp(model.type, 'bilbka')
        model.gait2d('Set','Bilateral BKA', model.anklestiffness, model.ankledamping);
	elseif strcmp(model.type, 'aka')
		model.gait2d('Set','Right AKA', model.anklestiffness, model.ankledamping);	
	end
	
	% ask what the mass is, and store it
	model.mass = model.gait2d('Get','Total Mass');
	% fprintf('Total mass of model: %8.3f kg\n',model.mass);
	
	% ask for the Fmax values of the muscles
	model.Fmax = model.gait2d('Get','Fmax');

end
