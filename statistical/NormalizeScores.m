function [Scores_norm]=NormalizeScores(allscores)
allscores_=allscores+max(abs(allscores));
Scores_norm=(allscores_)./max((allscores_));
end