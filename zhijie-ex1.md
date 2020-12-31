-- R = QQ[a_1,b_1,c_1,a_2,b_2,a_3,b_3,a_4,b_4,a_5,b_5,a_6,b_6,Degrees=>{{1,0,0,0},{1,0,0,1},{1,0,0,2},{1,1,0,0},{1,1,0,1},{1,1,1,0},{1,1,1,1},{0,1,0,0},{0,1,0,1},{0,1,1,0},{0,1,1,1},{0,0,1,0},{0,0,1,1}}]

R = QQ[a_1,b_1,c_1,a_2,b_2,a_3,b_3,a_4,b_4,a_5,b_5,a_6,b_6,Degrees=>{{1,0,0},{1,0,0},{1,0,0},{1,1,0},{1,1,0},{1,1,1},{1,1,1},{0,1,0},{0,1,0},{0,1,1},{0,1,1},{0,0,1},{0,0,1}}]

J = ideal(b_1,a_6,a_4,a_2,a_1,-a_5*b_2+a_3*b_4,b_2*b_6+a_3,b_4*b_6+a_5)

hJnu = hilbertFunction({1,1,1},J)

-- basis({1,1,1},R/J) --- figure out how to read this; output is matrix {{c_1*b_4*b_6, c_1*b_5, b_2*b_6, b_3}}
-- compare matrix entries with the list:
-- fL = c_1*a_5,a_3,c_1*b_4*b_6,b_2*b_6,c_1*a_5+c_1*b_4*b_6,c_1*b_5,a_3+c_1*b_5+b_2*b_6,b_3
-- which is got by applying f to l in L = apply(indsnu,i->p_i)
-- note that ideal(fL):J = ideal(c_1*b_4*b_6,b_2*b_6,c_1*b_5,c_1*a_5,b_3,a_3) which is almost the same as the matrix above

Jnu = ideal(c_1*a_5,a_3,c_1*b_4*b_6,b_2*b_6,c_1*a_5+c_1*b_4*b_6,c_1*b_5,a_3+c_1*b_5+b_2*b_6,b_3)

hilbertSeries(Jnu,Reduce=>true)

hilbertSeries(J,Reduce=>true)

plm = matrix{{1,0,-c_1,0,0,-b_2,0,0,0,-a_3,-b_3,0,0,0}, {0,1,0,-c_1,0,0,-b_2,0,0,0,-a_3,-b_3,0,0},{0,0,0,1,0,-b_4,0,0,0,-a_5,-b_5,0,0,0},{0,0,0,0,1,0,-b_4,0,0,0,-a_5,-b_5,0,0},{0,0,0,0,0,0,1,0,0,0,-b_6,0,0,0},{0,0,0,0,0,0,0,1,0,0,0,-b_6,0,0},{0,0,0,0,0,0,0,0,1,0,0,0,-b_6,0},{0,0,0,0,0,0,0,0,0,0,0,1,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,1,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,1}}

s = subsets(14,10);

e0 = {0..9}

mins := {}

for e in s do if ( det(plm^e0_e) != 0 ) then mins = append(mins,det(plm^e0_e))

inds := {}

for e in s do if ( det(plm^e0_e) != 0 ) then inds = append(inds,e)

wts := {} 

for e in s do if ( det(plm^e0_e) != 0 ) then wts = append(wts,degree(det(plm^e0_e)))

pluck = QQ[apply(inds,i->p_i)]

f = map(R,pluck,mins)

Q = R/J

fbar = map(Q,pluck,mins)

K = kernel fbar

Kh = homogenize(K,p_{0,1,3,4,6,7,8,11,12,13})

indsnu := {}

-- for i in 0..#wts-1 do if take(wts#i,3)=={1,1,1} then indsnu = append(indsnu,inds#i)
for i in 0..#wts-1 do if wts#i=={1,1,1} then indsnu = append(indsnu,inds#i)

I = ideal(apply(indsnu,i->p_i))

Knu = K + I

dim Knu

Knuh = homogenize(Knu,p_{0,1,3,4,6,7,8,11,12,13})

dim Knuh