classdef PGM_factor
    properties
        scope
        potentials
    end
    methods
        function obj = PGM_factor(scope, potentials)
            if nargin < 2
                obj = PGM_factor([],1);
                return;
            end
            obj.scope = scope; % requires scope to be in order
            obj.potentials = potentials;
        end
        function objC = times(objA,objB)
            sampleSpaceA = size(objA.potentials);
            sampleSpaceA(sampleSpaceA == 1) = [];
            sampleSpaceB = size(objB.potentials);
            sampleSpaceB(sampleSpaceB == 1) = [];
            combinedScope = union(objA.scope,objB.scope);
            newShapeA = ones(1,numel(combinedScope));
            
            membersOfA = ismember(combinedScope,objA.scope);
            
            
            newShapeA(ismember(combinedScope,objA.scope)) = sampleSpaceA;
            
            membersOfB = ismember(combinedScope,objB.scope);
            
            newShapeB = ones(1,numel(combinedScope));
            newShapeB(ismember(combinedScope,objB.scope)) = sampleSpaceB;
            
            repA = newShapeB;
            repA(membersOfA & membersOfB) = 1;
            repB = newShapeA;
            repB(membersOfA & membersOfB) = 1;
            
            if numel(newShapeA) == 1
                newShapeA = [1,newShapeA];
                repB = [1,repB];
            end
            if numel(newShapeB) == 1
                newShapeB  = [1,newShapeB];
                repA = [1,repA];
            end
            
            newA = repmat(reshape(objA.potentials,newShapeA),repA);
            newB = repmat(reshape(objB.potentials,newShapeB),repB);
            objC = PGM_factor(combinedScope,newA.*newB);
        end
        function objC = mtimes(objA,objB)
            objC = times(objA,objB);
        end
        function objC = plus(objA,objB)
            sampleSpaceA = size(objA.potentials);
            sampleSpaceA(sampleSpaceA == 1) = [];
            sampleSpaceB = size(objB.potentials);
            sampleSpaceB(sampleSpaceB == 1) = [];
            combinedScope = union(objA.scope,objB.scope);
            newShapeA = ones(1,numel(combinedScope));
            
            membersOfA = ismember(combinedScope,objA.scope);
            
            
            newShapeA(ismember(combinedScope,objA.scope)) = sampleSpaceA;
            
            membersOfB = ismember(combinedScope,objB.scope);
            
            newShapeB = ones(1,numel(combinedScope));
            newShapeB(ismember(combinedScope,objB.scope)) = sampleSpaceB;
            
            repA = newShapeB;
            repA(membersOfA & membersOfB) = 1;
            repB = newShapeA;
            repB(membersOfA & membersOfB) = 1;
            
            if numel(newShapeA) == 1
                newShapeA = [1,newShapeA];
            end
            if numel(newShapeB) == 1
                newShapeB  = [1,newShapeB];
            end
            
            newA = repmat(reshape(objA.potentials,newShapeA),repA);
            newB = repmat(reshape(objB.potentials,newShapeB),repB);
            objC = PGM_factor(combinedScope,newA+newB);
        end
        function objB = marginalize(objA,nodesToRemove)
            scopeB = objA.scope;
            potentialsB = objA.potentials;
            sampleSpace = size(potentialsB);
            sampleSpace(sampleSpace == 1) = [];
            for ii = nodesToRemove
                ind = find(scopeB == ii);
                potentialsB = sum(potentialsB,ind);
            end
            
            [scopeB,kept,~] = setxor(scopeB,nodesToRemove);
            if numel(kept) > 1
                potentialsB = reshape(potentialsB,sampleSpace(kept));
            else
                potentialsB = potentialsB(:)';
            end
            objB = PGM_factor(scopeB,potentialsB);                
        end
        function objB = substitution(objA,chosenNodes,values)
            [scopeB,kept,~] = setxor(objA.scope,chosenNodes);
            withinScope = ismember(scopeB,objA.scope);
            scopeB = scopeB(withinScope);
            sampleSpace = size(objA.potentials);
            sampleSpace(sampleSpace == 1) = [];
            inds = cell(numel(sampleSpace),1);
            for ii = 1:numel(inds)
                chosenOne = find(chosenNodes == objA.scope(ii));
                if isempty(chosenOne)
                    inds{ii} = 1:sampleSpace(ii);
                else
                    inds{ii} = values(chosenOne);
                end
            end
            if numel(inds) > 1
                potentialsB = cell2sub(objA.potentials,inds{:});
            elseif numel(inds) == 1
                potentialsB = objA.potentials(inds{1});
            else
                potentialsB = 1;
            end
            if numel(kept) > 1
                potentialsB = reshape(potentialsB,sampleSpace(kept));
            else
                potentialsB = potentialsB(:)';
            end
            objB = PGM_factor(scopeB,potentialsB);
        end
        function obj = normalize(obj)
            partitionFunction = sum(obj.potentials(:));
            if partitionFunction ~= 0
                obj.potentials = obj.potentials/partitionFunction;
            end
        end
    end
end
function output = cell2sub(X,varargin)
output = X(varargin{:});
end