function options = setdefaultoptions(options, defaults)
% Assign the defaults value to any field not specified in the strucy
% options
%
%  options = ManageOptions(defaults, options)
%
% INPUT
%  options (struct) the current options
%  defaults (struct) the default options
%
% OUTPUT
%  options (struct) the updated options
%

if nargin<2 || isempty(options)
    options = defaults;
else
    fields = fieldnames(options);
    defaultFields = fieldnames(defaults);
    numFields = length(defaultFields);
    for i=1:numFields
        if ~isfield(options, defaultFields(i))
            value = getfield(defaults, char(defaultFields(i)));
            options = setfield(options, char(defaultFields(i)), value);
        end
    end
end
