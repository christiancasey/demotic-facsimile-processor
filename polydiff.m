function diff = polydiff(p1, p2)

	pu = union(p1,p2);
	pi = intersect(p1,p2);
	
	au = area(pu);
	ai = area(pi);
	diff = (au- ai)/au;