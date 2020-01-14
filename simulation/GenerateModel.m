function [model, data] = GenerateModel(data)

switch data.nDoF
    case 9
        [model, data] = GenerateModel_9(data);
    case 12
        [model, data] = GenerateModel_12(data);
    case 42
        [model, data] = GenerateModel_42(data);
    otherwise
        disp(['No model for degree of freedom ', num2str(data.nDoF)])
end

end