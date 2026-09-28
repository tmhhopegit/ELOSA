classdef ExhaustiveRegionFinder   % original code; only this class name was fixed (it said pls_ensemble_args)
    properties
        
        Attempts = 0;
        NumTried = 0;
        NumFiles = 0;
        ParetoLength = 0;
        Groups = [];
        DiscreteGroups = [];
        DiscreteGroupSizes = [];
        DiscreteRegSizes = [];
        ParetoFront = [];
        NumLesions = 0;
        Lesions = [];
        GroupSizes = [];
        TargNum = 0;
        ChunkSize = 1000
        LastGroup = [];
    end
    methods
        function obj = ExhaustiveRegionFinder(Lesions)
            obj.Lesions = Lesions;
            obj.NumLesions = size(Lesions,1);
            %obj = obj.CalcGroupsByVoxel();
        end
        function obj = CalcGroupsByVoxel(obj)
            [obj.DiscreteGroups,~,IC] = unique(obj.Lesions','rows');
            obj.DiscreteRegSizes = zeros(size(obj.DiscreteGroups,1),1,'single');
%             for i=1:size(obj.DiscreteGroups,1)
%                 obj.DiscreteRegSizes(i) = length(find(IC==i));
%                 %disp([num2str(i) '/' num2str(size(obj.DiscreteGroups,1))])
%             end
            disp(['Found ' num2str(size(obj.DiscreteGroups,1))]);
            %obj.DiscreteGroupSizes = sum(single(obj.DiscreteGroups),2);
            obj.DiscreteGroups = logical(obj.DiscreteGroups);
        end
        function obj = GrowGroup(obj,GroupToGrow,LastItemAdded,MinRegSize,KeyGroups,TabuList)
            if(obj.NumTried>=obj.ChunkSize)
                obj = obj.RecordGrowing(GroupToGrow,MinRegSize);
            end
            obj.NumTried = obj.NumTried + 1;
            obj.Groups(obj.NumTried,:) = GroupToGrow;
            skg = size(KeyGroups,1);
            s = sum(KeyGroups(:,(LastItemAdded+1):obj.NumLesions),1);
            s = find(s>1 & s<skg) + LastItemAdded;
            for i=s
                KeyGp = KeyGroups(KeyGroups(:,i),:);
                TestGp = GroupToGrow; TestGp(i) = true;
                AllIn = all(KeyGp,1);
                if(isempty(find(AllIn(~TestGp(1:i)),1)))
                    obj = obj.GrowGroup(AllIn,i,MinRegSize,KeyGp);
                end
            end
        end
        function obj = GrowGroup_v2(obj,GroupToGrow,LastItemAdded,MinRegSize,KeyGroups,Remaining)
            if(obj.NumTried>=obj.ChunkSize)
                obj = obj.RecordGrowing(GroupToGrow,MinRegSize);
            end
            obj.NumTried = obj.NumTried + 1;
            obj.Groups(obj.NumTried,:) = GroupToGrow;
            After = Remaining > LastItemAdded;
            for i=find(After,1):length(Remaining)
                KeyGp = KeyGroups(KeyGroups(:,i),:);
                TestGp = GroupToGrow; TestGp(Remaining(i)) = true;
                AllIn = all(KeyGp,1);
                if(isempty(find(AllIn(1:(i-1)),1)))
                    sums = sum(KeyGp,1);
                    TestGp(Remaining) = AllIn;
                    KeyR = Remaining;
                    Remove = AllIn | sums<=1;
                    KeyR(Remove) = [];
                    KeyGp(:,Remove) = [];
                    obj = obj.GrowGroup_v2(TestGp,Remaining(i),MinRegSize,KeyGp,KeyR);
                end
            end
        end
        function obj = GrowGroup_v3(obj,GroupToGrow,LastItemAdded,MinRegSize,KeyGroups,Remaining)
            if(obj.NumTried>=obj.ChunkSize)
                obj = obj.RecordGrowing(GroupToGrow,MinRegSize);
            end
            obj.NumTried = obj.NumTried + 1;
            obj.Groups(obj.NumTried,:) = GroupToGrow;
            [uKG,x,y] = unique(KeyGroups','rows'); uKG=uKG';
            After = Remaining > LastItemAdded;
            y = [y,(1:length(y))'];
            KeyCombs = y(After,:); KeyCombs = KeyCombs(~ismember(KeyCombs(:,1),y(~After,1)),:);
            [~,b]=unique(KeyCombs(:,1),'stable');
            KeyCombs = KeyCombs(b,:);
            for i=1:size(KeyCombs,1)
                KeyGp = KeyGroups(uKG(:,KeyCombs(i,1)),:);
                AllIn = all(KeyGp,1);
                if(isempty(find(AllIn(1:(KeyCombs(i,2)-1)),1)))
                    sums = sum(KeyGp,1);
                    TestGp = GroupToGrow; 
                    TestGp(Remaining) = AllIn;
                    KeyR = Remaining;
                    Remove = AllIn | sums<=1;
                    KeyR(Remove) = [];
                    KeyGp(:,Remove) = [];
                    obj = obj.GrowGroup_v3(TestGp,Remaining(KeyCombs(i,2)),MinRegSize,KeyGp,KeyR);
                end
            end
        end
        function obj = GrowGroup2(obj,GroupToGrow,LastItemAdded,MinRegSize,KeyGroups,Inds,TabooList)
            if(obj.NumTried >= obj.ChunkSize)
                obj.NumFiles = obj.NumFiles + 1;
                fname = ['L14_Interim',num2str(MinRegSize),'_',num2str(obj.NumFiles),'.dat'];
                x = sparse(obj.Groups);
                save(fname,'x');
                x = [];
                g2 = find(obj.Groups(size(obj.Groups,1),:));
                obj.Groups = false(obj.ChunkSize,obj.NumLesions);
                obj.NumTried = 1;
                Line = '';
                g = find(GroupToGrow);
                for i=1:length(g)
                    Line = [Line num2str(g(i)) ','];
                end
                Line(length(Line)) = [];
                x = toc();
                disp(['P: ' num2str(MinRegSize) '; N: ' num2str(obj.Attempts) '; T: ' num2str(x) '; G: ' Line]);
                obj.Attempts = obj.Attempts + obj.ChunkSize;
                g1 = find(obj.LastGroup); 
                Found = 0; ind=1;
                while(Found==0 && ind < length(g1) && ind < length(g2))
                    if(g1(ind)>g2(ind))
                        Found = 1;
                    elseif(g1(ind)<g2(ind))
                        Found = -1;
                    end
                    ind = ind + 1;
                end
                obj.LastGroup = GroupToGrow;
                tic;
            end
            [a,~,c]=unique(KeyGroups','rows');
            a = a';
            for i=1:size(a,2)
                KeyGp = KeyGroups(a(:,i),:);
                Gp = all(KeyGp,1);
                Gpi = Inds(Gp);
                if(Gpi(1)>LastItemAdded)
                    NewGroup = GroupToGrow;
                    NewGroup(Gpi) = true;
                    obj.NumTried = obj.NumTried + 1;
                    obj.Groups(obj.NumTried,:) = NewGroup;
                    NoLesion = sum(KeyGp,1) <= 1;
                    tempInds = Inds;
                    Irrelevant = (NoLesion | Gp);
                    %tempInds(Irrelevant) = [];
                    if(size(KeyGp,1)>2)
                        %KeyGp = KeyGp(:,~Irrelevant);
                        NewLast = find(c==i,1);
                        obj = obj.GrowGroup2(NewGroup,NewLast,MinRegSize,KeyGp,tempInds);
                    end
                end
            end
        end
        function obj = RecordGrowing(obj,GroupToGrow,MinRegSize)
            try
                obj.NumFiles = obj.NumFiles + 1;
                fname = ['L18_Interim',num2str(MinRegSize),'_',num2str(obj.NumFiles),'.dat'];
                x = sparse(obj.Groups);
                save(fname,'x');
                x = [];
                if(false && ~isempty(obj.LastGroup))
                    g = find(GroupToGrow);
                    gl = find(obj.LastGroup);
                    FirstDiff = 1;
                    while(g(FirstDiff)==gl(FirstDiff))
                        FirstDiff = FirstDiff + 1;
                    end
                else
                    g = find(GroupToGrow);
                    FirstDiff = length(g);
                end
                Line = '';
                for i=1:FirstDiff
                    Line = [Line num2str(g(i)) ','];
                end
                Line(length(Line)) = [];
                if(~isempty(obj.LastGroup) && all(obj.LastGroup==GroupToGrow))
                    disp('Error!')
                end
                obj.LastGroup = GroupToGrow;
                
                x = toc();
                obj.Attempts = obj.Attempts + obj.ChunkSize;
                disp(['P: ' num2str(MinRegSize) '; L: ' num2str(MinRegSize) '; T: ' num2str(x) '; G: ' Line]);
                obj.Groups = false(obj.ChunkSize,obj.NumLesions);
                obj.NumTried = 0;
                tic;
            catch ME
                disp(ME.message)
            end
        end
        function obj = GrowGroup3(obj,GroupToGrow,LastItemAdded,MinRegSize,KeyGroups,RemainingInds)
            if(obj.NumTried >= obj.ChunkSize)
                obj = obj.RecordGrowing(GroupToGrow,MinRegSize);
            end
            s = size(KeyGroups,1);
            obj.NumTried = obj.NumTried + 1;
            obj.Groups(obj.NumTried,:) = GroupToGrow;
            if(s>2)
                a = unique(KeyGroups','rows'); a = a';
                for i=1:size(a,2)
                    KeyGp = KeyGroups(a(:,i),:);
                    AllIn = all(KeyGp,1);
                    IndsI = RemainingInds(AllIn);
                    if(IndsI(1) > LastItemAdded)
                        NewGroup = GroupToGrow; NewGroup(IndsI) = true;
                        NoLesion = all(~KeyGp,1);
                        KeyGp(:,NoLesion | AllIn) = [];
                        KeyInds = RemainingInds; KeyInds(NoLesion | AllIn) = [];
                        obj = obj.GrowGroup3(NewGroup,IndsI(1),MinRegSize,KeyGp,KeyInds);
                    end
                end
            end
        end
        function obj = GrowGroup4(obj,GroupToGrow,LastItemAdded,MinRegSize,KeyGroups,RemainingInds)
            if(obj.NumTried >= obj.ChunkSize)
                obj = obj.RecordGrowing(GroupToGrow,LastItemAdded);
            end
            sk = size(KeyGroups,1);
            obj.NumTried = obj.NumTried + 1;
            obj.Groups(obj.NumTried,:) = GroupToGrow;
            if(sk>2)
                s=sum(KeyGroups,1);
                pot = find(s>1 & s<sk);
%                 [a,b,c] = unique(KeyGroups','rows');
%                 a = a';
%                 x = sum(a,1)>1;
%                 a = a(:,x);
%                for i=1:size(a,2)
                for i=1:length(pot)
                    KeyGp = KeyGroups(KeyGroups(:,pot(i)),:);
                    sums = sum(KeyGp,1);
                    AllIn = sums==size(KeyGp,1);
                    IndsI = RemainingInds(AllIn);
                    if(IndsI(1) > LastItemAdded)
                        NewGroup = GroupToGrow; NewGroup(IndsI) = true;
                        NoLesion = sums<=1; KeyInds = RemainingInds;
%                         KeyGp(:,NoLesion | AllIn) = [];
%                         KeyInds = RemainingInds; KeyInds(NoLesion | AllIn) = [];
                        obj = obj.GrowGroup4(NewGroup,pot(i),MinRegSize,KeyGp,KeyInds);
                    end
                end
            end
        end
        function obj = FindRecursive(obj,MinRegSize,Indices)
            obj.ParetoLength = 0;
            obj.NumTried = 0;
            obj.ChunkSize = 1000000;
            %ParetoGraph = nan(obj.NumLesions,1);
            obj.Groups = false(obj.ChunkSize,obj.NumLesions);
            obj.ParetoFront = false(obj.ChunkSize,1);
            obj.Attempts = 0;
            if(nargin < 3 || isempty(Indices))
                Indices = (1:obj.NumLesions)';
            end
            tic;
            for i=1:length(Indices)
                p = Indices(i);
                obj.Attempts = obj.Attempts + 1;
                Gp = false(1,obj.NumLesions);
                Gp(p) = true;
                %Gp = sparse(Gp);
                KeyGroups = (obj.DiscreteGroups(all(obj.DiscreteGroups(:,Gp),2),:));
                if(size(KeyGroups,1)>1)
                    Gp = all(KeyGroups,1);
                    fGp = find(Gp);
                    if(fGp(1)>=Indices(i))
                        KeyG = KeyGroups; KeyG(:,Gp) = []; KeyR = 1:obj.NumLesions; KeyR(Gp) = [];
                        obj = obj.GrowGroup_v3(Gp,p,MinRegSize,KeyG,KeyR);
                    end
                end
%                 x = obj.Groups (obj.ParetoFront,:);
%                 n = ['BigGroups_',num2str(Indices(1)),'.dat'];
%                 dlmwrite(n,x);
%                 x = [];
%                 ParetoGraph(p,1) = length(find(obj.ParetoFront));
%                 plot(ParetoGraph);
%                 drawnow;
%                 pause(0.1)
            end
            obj.Groups = obj.Groups(1:obj.NumTried,:);
            obj.Groups = unique(obj.Groups,'rows');
        end
        function obj = FindRecursive2(obj,MinRegSize,Indices)
            obj.ParetoLength = 0;
            obj.NumTried = 0;
            obj.ChunkSize = 10000;
            %ParetoGraph = nan(obj.NumLesions,1);
            obj.Groups = false(obj.ChunkSize,obj.NumLesions);
            obj.ParetoFront = false(obj.ChunkSize,1);
            obj.Attempts = 0;
            if(nargin < 3 || isempty(Indices))
                Indices = (1:obj.NumLesions)';
            end
            tic;
            for i=1:length(Indices)
                p = Indices(i);
                obj.Attempts = obj.Attempts + 1;
                Gp = false(1,obj.NumLesions);
                Gp(p) = true;
                %Gp = sparse(Gp);
                KeyGroups = (obj.DiscreteGroups(all(obj.DiscreteGroups(:,Gp),2),:));
                if(size(KeyGroups,1)>1)
                    Gp = all(KeyGroups,1);
                    fGp = find(Gp);
                    if(fGp(1)>=Indices(i))
                        Inds = 1:size(KeyGroups,2);
                        NoLesion = sum(KeyGroups,1)<=1;
                        Inds(Gp|NoLesion) = [];
                        KeyGroups = KeyGroups(:,Inds);
%                         obj.NumTried = obj.NumTried+1;
%                         obj.Groups(obj.NumTried,:) = Gp;
                        obj = obj.GrowGroup4(Gp,p,MinRegSize,KeyGroups,Inds);
                    end
                end
%                 x = obj.Groups (obj.ParetoFront,:);
%                 n = ['BigGroups_',num2str(Indices(1)),'.dat'];
%                 dlmwrite(n,x);
%                 x = [];
%                 ParetoGraph(p,1) = length(find(obj.ParetoFront));
%                 plot(ParetoGraph);
%                 drawnow;
%                 pause(0.1)
            end
            obj.Groups = obj.Groups(1:obj.NumTried,:);
            obj.Groups = unique(obj.Groups,'rows');
        end
        function obj = FindRecursiveSeeded(obj,Seeds)
            obj.ParetoLength = 0;
            obj.NumTried = 0;
            obj.ChunkSize = 10000;
            %ParetoGraph = nan(obj.NumLesions,1);
            obj.Groups = false(obj.ChunkSize,obj.NumLesions);
            obj.ParetoFront = false(obj.ChunkSize,1);
            obj.Attempts = obj.ChunkSize;
            Last = Seeds(:,1);
            Seeds(:,1) = [];
            Seeds = logical(Seeds);
            tic;
            for i=1:size(Seeds,1)
                obj.Attempts = obj.Attempts + 1;
                Gp = false(1,obj.NumLesions);
                Gp(Seeds(i,:)) = true;
                p = Last;
                KeyGroups = (obj.DiscreteGroups(all(obj.DiscreteGroups(:,Gp),2),:));
                if(size(KeyGroups,1)>1)
                    Gp = all(KeyGroups,1);
                    obj.NumTried = obj.NumTried+1;
                    obj = obj.GrowGroup(Gp,p,i*-1,KeyGroups,[]);
                end
            end
            obj.Groups = obj.Groups(sum(obj.Groups,2)>0,:);
            obj.Groups = unique(obj.Groups,'rows');
        end
        function obj = FindIncremental(obj)
            Chunk = 1000;
            gp = [obj.Lesions(1,:) ; ~obj.Lesions(1,:)];
            fnum = 1;
            save(['IncrementalGp_',num2str(fnum),'.dat'],'gp');
            NumGroups = 2;
            obj.Groups = logical(obj.Groups);
            obj.Lesions = logical(obj.Lesions);
            NumGp = 2;
            fnames = {['IncrementalGp_',num2str(fnum),'.dat']};
            for i=2:obj.NumLesions
                if(length(find(obj.Lesions(i,:)))>0)
                    if(fnum > 1)
                        NumF = fnum - 1;
                        fout = cell(NumF,1);
                        MoreNums = cell(NumF,1);
                        parfor j=1:NumF
                            gp = importdata(fnames{j});
                            KeyGp = gp(:,obj.Lesions(i,:));
                            keygp = find(~all(KeyGp,2) & ~all(~KeyGp,2));
                            newbit = gp(keygp,:) & repmat(obj.Lesions(i,:),length(keygp),1);
                            KeyGp = gp(:,~obj.Lesions(i,:));
                            keygp = find(~all(KeyGp,2) & ~all(~KeyGp,2));
                            newbit = [newbit ; gp(keygp,:) & repmat(~obj.Lesions(i,:),length(keygp),1)];
                            fout(j) = {['IncrementalGp_',num2str(fnum),'.dat']};
                            parsave(fout{j},newbit);
                            MoreNums(j) = {size(newbit,1)};
                            disp(num2str([i,j,size(newbit,1)]))
                        end
                        fnames = [fnames ; fout];
                        fnums = length(fnames);
                        NumGp = NumGp + sum(cell2mat(MoreNums));
                        disp(num2str([i,NumGp]))
                    else
                        KeyGp = gp(:,obj.Lesions(i,:));
                        keygp = find(~all(KeyGp,2) & ~all(~KeyGp,2));
                        newbit = gp(keygp,:) & repmat(obj.Lesions(i,:),length(keygp),1);
                        NumGp = NumGp + size(newbit,1);
                        gp = [gp ; newbit];
                        KeyGp = gp(:,~obj.Lesions(i,:));
                        keygp = find(~all(KeyGp,2) & ~all(~KeyGp,2));
                        newbit = gp(keygp,:) & repmat(obj.Lesions(i,:),length(keygp),1);
                        NumGp = NumGp + size(newbit,1);
                        gp = [gp ; newbit];
                        disp(num2str([i,NumGp]))
                    end
                end
                if(size(gp,1)>10000)
                    fnames(fnum) = {['IncrementalGp_',num2str(fnum),'.dat']};
                    save(fnames{fnum},'gp');
                    fnum = fnum+1;
                    gp = [];
                end
            end
                
%             
%             obj.ParetoLength = 0;
%             ParetoGraph = [];
%             obj.Groups = false(obj.NumLesions,obj.NumLesions);
%             for i=1:obj.NumLesions
%                 obj.Groups(i,i) = true;
%             end
%             obj.ParetoFront = true(obj.NumLesions,1);
%             obj.Attempts = 0; ind = 1;
%             for i=1:obj.NumLesions
%                 
%                 tic;
%                 ParetoGraph(ind,1) = size(obj.Groups,1);
%                 ind = ind + 1;
%                 plot(ParetoGraph);
%                 drawnow;
%                 pause(0.1)
%                 
%                 TempGroups = obj.Groups;
%                 TempGroups(:,i) = true;
%                 TempGroups = TempGroups(TempGroups(:,i)==0,:);
%                 Possibles = obj.DiscreteGroups(obj.DiscreteGroups(:,i)>0,:);
%                 NotPossibles = sum(Possibles,1) == 0;
%                 NotPossibles = sum(TempGroups(:,NotPossibles),2) > 0;
%                 TempGroups = TempGroups(~NotPossibles,:);
%                  
%                 if(~isempty(TempGroups))
% 
%                     TempGroups(:,i) = true;
%                     d = obj.DiscreteGroups;
%                     DatCell = cell(16,1);
%                     GpCell = cell(16,1);
%                     for j=1:16
%                         inds = mod(1:size(TempGroups,1),16) + 1;
%                         DatCell(j,1) = {d};
%                         GpCell(j,1) = {TempGroups(inds==j,:)};
%                     end
%                     
%                     parfor j = 1:16
%                         d = DatCell{j,1};
%                         g = GpCell{j,1};
%                         DatCell(j,1) = {[]};
%                         GpCell(j,1) = {[]};
%                         
% %                         g = CheckNewGroups_mex(d,g);
% %                         GpCell(j,1) = {g};
%                         
%                         Exists = false(size(g,1),1);
%                         
%                         for k = 1:size(g,1)
%                             s = all(d(:,g(k,:)),2);
%                             g(k,:) = all(d(s,:),1);
%                             if(~isempty(find(s,1)))
%                                 Exists(k,1) = true;
%                                 g(k,:) = all(d(s,:),1);
%                             end
%                         end
%                         
%                         GpCell(j,1) = {g(Exists,:)};
%                     end
%                     InitLen = size(TempGroups,1);
%                     FinalLen = size(cell2mat(GpCell),1);
%                     TempGroups = unique([obj.Groups ; cell2mat(GpCell)],'rows');
%                     obj.Groups = TempGroups;
%                     disp(num2str([i,size(obj.Groups,1),InitLen,FinalLen,toc]))
%                     save 'temp.dat' 'TempGroups'
%                     
% 
%                 end
%             end
        end
        function obj = FindIncrementalSparse(obj)
            obj.ParetoLength = 0;
            ParetoGraph = [];
            obj.Groups = false(obj.NumLesions,obj.NumLesions);
            for i=1:obj.NumLesions
                obj.Groups(i,i) = true;
            end
            MainGroups = sparse(obj.Groups);
            DiscGroups = sparse(obj.DiscreteGroups);
            obj.ParetoFront = true(obj.NumLesions,1);
            obj.Attempts = 0; ind = 1;
            for i=1:obj.NumLesions
                
                tic;
                ParetoGraph(ind,1) = size(MainGroups,1);
                ind = ind + 1;
                plot(ParetoGraph);
                drawnow;
                pause(0.1)
                
                TempGroups = MainGroups;
                TempGroups = TempGroups(TempGroups(:,i)==0,:);
                Possibles = DiscGroups(DiscGroups(:,i)>0,:);
                NotPossibles = sum(Possibles,1) == 0;
                NotPossibles = sum(TempGroups(:,NotPossibles),2) > 0;
                TempGroups = TempGroups(~NotPossibles,:);
                 
                if(~isempty(TempGroups))

                    TempGroups(:,i) = true;
                    d = DiscGroups;
                    DatCell = cell(16,1);
                    GpCell = cell(16,1);
                    for j=1:16
                        inds = mod(1:size(TempGroups,1),16) + 1;
                        DatCell(j,1) = {d};
                        GpCell(j,1) = {TempGroups(inds==j,:)};
                    end
                    
                    parfor j = 1:16
                        d = DatCell{j,1};
                        g = GpCell{j,1};
                        DatCell(j,1) = {[]};
                        GpCell(j,1) = {[]};
                        
%                         g = CheckNewGroups_mex(d,g);
%                         GpCell(j,1) = {g};
                        
                        Exists = false(size(g,1),1);
                        
                        for k = 1:size(g,1)
                            s = all(d(:,g(k,:)),2);
                            if(~isempty(find(s,1)))
                                Exists(k,1) = true;
                                g(k,:) = all(d(s,:),1);
                            end
                        end
                        
                        GpCell(j,1) = {g(Exists,:)};
                    end
                    disp(num2str([i,size(MainGroups,1),size(TempGroups,1),toc]))
                    MainGroups = unique([MainGroups ; cell2mat(GpCell)],'rows');
                    save 'temp.dat' 'MainGroups'
                    

                end
            end
            obj.Groups = MainGroups;
        end
    end
end
