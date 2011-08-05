function [positives]=CountPositive(SQ)
found=regexp(SQ,'[RK]');
positives=length(found);

end