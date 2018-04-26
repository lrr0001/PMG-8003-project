classdef Gibbs_sampling
    properties
        factors
        blocks
        blockFactorCorr
        scope
        states
    end
    methods
        function obj = Gibbs_sampling(factors,blocks)
            if nargin < 2
                blocks = num2cell(multiUnion(factors(:).scope));
                obj = Gibbs_sampling(factors,blocks);
            end
            obj.factors = factors;
            obj.blocks = blocks;
            obj.scope = multiUnion(factors(:).scope);
            obj.blockFactorCorr = findBlockFactorCorr(factors,blocks);
            obj.states = cell(1,1);
            obj.states{1} = ones(numel(obj.scope),1);
        end
        function objA = updateStep(objA,blkNum,stateNum)
            if nargin < 3
                stateNum = 1;
            end
            [fixedVars,kept,~] = setxor(objA.scope,objA.blocks{blkNum});
            updateVars = setxor(objA.scope,objA.scope(kept));
            if numel(updateVars) == 0
                error('no variables to update');
            end
            fixedVals = objA.states{stateNum}(kept);
            transitionProb = PGM_factor([],1);
            for ii = objA.blockFactorCorr{blkNum}
                transitionProb = transitionProb*objA.factors(ii).substitution(fixedVars,fixedVals);
            end
            transitionProb = transitionProb.normalize();
            if numel(updateVars) ~= numel(transitionProb.scope)
                error('not all update variables appear in factors');
            end
            
            randomNum = rand;
            ind = find(randomNum >= 1 - cumsum(transitionProb.potentials(:)),1);
            newVals = cell(1,numel(size(transitionProb.potentials)));
            if numel(updateVars) > 1
                [newVals{:}] = ind2sub(size(transitionProb.potentials),ind);
                objA.states{stateNum}(updateVars) = cell2mat(newVals);
            elseif numel(updateVars) == 1
                objA.states{stateNum}(updateVars) = ind;
            end 
        end
        function obj = updateSweep(obj,stateNum)
            if nargin == 1
                stateNum = 1;
            end
            for ii = 1:numel(obj.blocks)
                obj = obj.updateStep(ii,stateNum);
            end
                
        end
        function samples = getGSamples(obj,nos,burnIn,corrGap,init)
            if nargin < 5
                init = obj.states{1};
            end
            obj.states{1} = init;
            for ii = 1:burnIn
                if mod(ii,40)
                    disp(sprintf('burn-in, index = %d',ii));
                end
                obj = obj.updateSweep(1);
            end
            samples = ones(numel(obj.scope),nos);
            for ii = 1:nos
                if mod(ii,20) == 0
                    ii
                end
                samples(:,ii) = obj.states{1};
                for jj = 1:corrGap
                    obj = obj.updateSweep(1);
                end
            end
            
        end
            
    end
end
function blockFactorCorr = findBlockFactorCorr(factors,blocks)
blockFactorCorr = cell(numel(blocks),1);
for ii = 1:numel(blocks)
    temp = [];
    for jj = 1:numel(factors)
        if any(ismember(factors(jj).scope,blocks{ii}))
            temp = [temp,jj];
        end
    end
    blockFactorCorr{ii} = temp;
end
end
function A = multiUnion(varargin)
A = [];
for ii = 1:nargin
    A = union(A,varargin{ii});
end
end