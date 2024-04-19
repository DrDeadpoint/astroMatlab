function ret = checkClass(in,desClass)
if ~isa(in,desClass)
    error(['Incorrect class: must be ' desClass])
end
end