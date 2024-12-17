import pickle
with open('3protpred.p', 'rb') as f:
    scores = pickle.load(f)
baseline = scores
l = list(scores.items())
def all_bigger(l1, l2, num):
    cnt = 0
    for i in range(5):
        if l1[i] > l2[i]:
            cnt += 1
    return cnt == num

def med_bigger(l1, l2, num):
    l1.sort()
    l2.sort()
    return sum(l1[1:4]) > sum(l2[1:4])

f = open('report.txt', '+w')

f.write('baseline:\n')
f.write('GRDLYDD: ' + str(baseline['GRDLYDD']) + '\n')

def iterate(cond, basel, num):
    for p in l:
        if cond(p[1], basel, num):
            f.write(p[0]+ ': ' + str(p[1]) + '\n')

f.write('all bigger:\n')
iterate(all_bigger, baseline['GRDLYDD'], 5)
f.write('4 bigger:\n')
iterate(all_bigger, baseline['GRDLYDD'], 4)
f.write('med bigger:\n')
iterate(med_bigger, baseline['GRDLYDD'], 0)
f.close()
