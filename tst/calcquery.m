function ans = calcquery(fn,dt)
	fid = fopen(fn,'rt');
	p0 = readvec(fid);
	Q = readmat(fid);
	ns = size(Q,1);
	[times,types,inds] = readtraj(fid);
	tpts = times(1):dt:times(end);
	A = expm(Q*dt);
	alpha = genalpha(p0,Q,A,times,types,inds,tpts);
	beta = genbeta(p0,Q,A,times,types,inds,tpts);
	smooth = alpha.*beta;
	smoothsum = sum(smooth,1);
	smooth = smooth./repmat(smoothsum,[ns 1]);
	ans = [];
	while(1)
		qtype = fscanf(fid,'%d',1);
		if (qtype==-1 || feof(fid)) break;
		end;
		if (qtype==0)
			inds1 = readvec(fid);
			inds1 = inds1+1;
			t = fscanf(fid,'%f',1);
			a = findclose(compress(alpha,inds1),tpts,t);
			ans = [ans;a];
		elseif (qtype==1)
			inds1 = readvec(fid);
			inds1 = inds1+1;
			t = fscanf(fid,'%f',1);
			a = findclose(compress(smooth,inds1),tpts,t);
			ans = [ans;a];
		elseif (qtype==2)
			inds1 = readvec(fid);
			inds1 = inds1+1;
			a = sum(compress(smooth,inds1))...
					/length(tpts)*(times(end)-times(1));
			ans = [ans;a];
		elseif (qtype==3)
			inds1 = readvec(fid);
			inds2 = readvec(fid);
			inds1 = inds1+1;
			inds2 = inds2+1;
			bp = permute(beta(inds2,:,:),[3 1 2]);
			ap = alpha(inds1,:)./repmat(smoothsum,[size(inds1,1) 1]);
			ss = (ap'*Q(inds1,inds2))';
			ss = ss.*beta(inds2,:);
			a = sum(sum(ss))...
					/length(tpts)*(times(end)-times(1));
			for i=2:length(times) % assume which variable is known
				if (types(i) && ...
						any(ismember(inds{i-1},inds1)) && ...
						~any(ismember(inds{i-1},inds2)) && ...
						any(ismember(inds{i},inds2)) && ...
						~any(ismember(inds{i},inds1)))
					a = a+1;
				end;
			end;
			ans = [ans;a];
		else
			return;
		end;
	end;
end;

function v = compress(vs,inds)
	v = sum(vs(inds,:),1);
end;

function x = findclose(vs,ts,t)
	[delta,i] = min(abs(ts-t));
	x = vs(i);
end;

function alpha = genalpha(p0,Q,A,times,types,inds,tpts)
	ns = size(Q,1);
	nt = length(tpts);
	alpha = zeros(ns,nt);
	v = zeros(ns,1);
	v(inds{1}) = p0(inds{1});
	v = v/sum(v);
	t = times(1);
	alpha(:,1) = v;
	ai = 1;
	ti = 1;
	while(ti<length(times) && ai<length(tpts))
		newv = A'*v;
		v = zeros(ns,1);
		v(inds{ti}) = newv(inds{ti});
		while (ti<length(times) && times(ti+1)<=tpts(ai+1))
			if (types(ti+1)) % var change
				cQ = zeros(ns,ns);
				cQ(inds{ti},inds{ti+1}) = Q(inds{ti},inds{ti+1});
				cQ = cQ-diag(diag(cQ));
				v = cQ'*v;
			end;
			ti = ti+1;
		end;
		v = v/sum(v);
		ai = ai+1;
		alpha(:,ai) = v;
	end;
end;

function beta= genbeta(p0,Q,A,times,types,inds,tpts)
	ns = size(Q,1);
	nt = length(tpts);
	beta = zeros(ns,nt);
	v = ones(ns,1);
	v = v/sum(v);
	t = times(1);
	bi = nt;
	beta(:,bi) = v;
	ti = length(times)-1;
	while(ti>0 && bi>1)
		while (ti>0 && times(ti)>tpts(bi-1))
			if (types(ti)) % var change
				cQ = zeros(ns,ns);
				cQ(inds{ti-1},inds{ti}) = Q(inds{ti-1},inds{ti});
				cQ = cQ-diag(diag(cQ));
				v = cQ*v;
			end;
			ti = ti-1;
		end;
		if (ti==0)
			ti = 1;
		end;
		newv = A*v;
		v = zeros(ns,1);
		v(inds{ti}) = newv(inds{ti});
		v = v/sum(v);
		bi = bi-1;
		beta(:,bi) = v;
	end;
end;

function [times,types,inds] = readtraj(fid)
	n = fscanf(fid,'%d',1);
	times = zeros(n+1,1);
	types = zeros(n+1,1);
	inds = cell(n,1);
	for i=1:n
		times(i) = fscanf(fid,'%f',1);
		types(i) = fscanf(fid,'%d',1);
		m = fscanf(fid,'%d',1);
		inds{i} = fscanf(fid,'%d',m);
		inds{i} = inds{i}+1;
	end;
	times(n+1) = fscanf(fid,'%f',1);
end;
function v = readvec(fid)
	n = fscanf(fid,'%f',1);
	v = fscanf(fid,'%f',n);
end;
function M = readmat(fid)
	m = fscanf(fid,'%f',1);
	n = fscanf(fid,'%f',1);
	M = fscanf(fid,'%f',m*n);
	M = reshape(M,m,n);
	M = M';
end;
	
