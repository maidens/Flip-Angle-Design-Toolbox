function check_system_requirements(  )
    %FUNCTION CHECK_SYSTEM_REQUIREMENTS verifies that the required MATLAB toolboxes are installed 
    % 
    %   Flip Angle Design Toolbox 
    %   John Maidens (maidens@eecs.berkeley.edu)
    %   June 2014

    v = ver;
    toolboxes = setdiff({v.Name}, 'MATLAB'); 
    if ~strncmp('Optimization Toolbox', toolboxes, 20)
        error('Flip Angle Design Toolbox requires MATLAB Optimization Toolbox')
    elseif ~strncmp('Symbolic Math Toolbox', toolboxes, 21)
        error('Flip Angle Design Toolbox requires MATLAB Symbolic Math Toolbox')
    end
    
end

