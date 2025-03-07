function res = structure_variables(structure, fieldname, variablename, value)
  if ~isfield(structure, fieldname)
    % If it doesn't exist, initialize it as an empty cell array
    structure.(fieldname).(variablename) = {};
    structure.(fieldname).(variablename){1} = value ;
  elseif ~isfield(structure.(fieldname), variablename)
    structure.(fieldname).(variablename) = {};
    structure.(fieldname).(variablename){1} = value ;
  else 
      structure.(fieldname).(variablename){end+1} = value ;
  end
  res = structure ;