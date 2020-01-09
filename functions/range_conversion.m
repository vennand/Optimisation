function new_value = range_conversion(old_value, old_max, old_min, new_max, new_min)
old_range = (old_max - old_min);
new_range = (new_max - new_min);
new_value = (((old_value - old_min) * new_range) / old_range) + new_min;
end