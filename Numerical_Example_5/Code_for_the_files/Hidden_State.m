function [State_Space] = Hidden_State(varargin)
%varargin should be a list of row vectors indicating the range of variation
%for each species
State_Space=allcomb(varargin{:});
end

