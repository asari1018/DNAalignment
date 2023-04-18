"""
Created on Thu Jul  9 18:41:45 2020

@author: arisa
"""
#2つの文字についてのスコアを計算する
def pair_score(str_1, str_2):
    if str_1 != '-' and str_2 != '-':
        if str_1 == str_2:
            score = 2
        else:
            score = -1
    else:
        if str_1 == '-' and str_2 == '-':
            score = 0
        else:
            score = -1
    
    #2つの文字についてのスコアを返り値とする
    return score

#二つの遺伝子配列のアラインメントのスコアを計算する
def calculate(seq_1, seq_2):
    #DP表を作成する
    scores = [[0 for j in range(len(seq_2)+1)] 
                for i in range(len(seq_1)+1)]
    #DP表を1マスずつ埋めていく
    for i in range(1,len(seq_1)+1):
        for j in range(1,len(seq_2)+1):
            #今回は、同じ文字同士->2点
            #異なる文字同士->-1点
            #GAP->-1点として、DP表の各座標が最大値をとるように
            #漸化式をたてて調べる
            tmp_score = pair_score(seq_1[i-1], seq_2[j-1])
            scores[i][j] = max([scores[i-1][j]-1, 
                                scores[i][j-1]-1, 
                                scores[i-1][j-1]+tmp_score])
    #二つの遺伝子配列のアラインメントのスコアを返り値とする
    return scores[len(seq_1)][len(seq_2)]

#二つの遺伝子配列をアラインメントする関数
def alignment(main_seq, sub_seq):
    #アラインメント後の配列を保存する配列を用意
    changed_sub = []
    changed_main = []
    #スコアのDP表を初期化
    scores = [[0 for j in range(len(sub_seq)+1)] 
                for i in range(len(main_seq)+1)]
    #それぞれのスコアがどこから計算されたかを記憶する配列を初期化
    how = [[0 for j in range(len(sub_seq)+1)] 
                for i in range(len(main_seq)+1)]

    #DP表の各座標が最大値をとるように漸化式をたてて調べる
    for i in range(1,len(main_seq)+1):
        for j in range(1,len(sub_seq)+1):
            tmp_score = pair_score(main_seq[i-1], sub_seq[j-1])
            scores[i][j] = max([scores[i-1][j]-1, 
                                scores[i][j-1]-1,
                                scores[i-1][j-1]+tmp_score])

            #ここで、後で配列を出力するため
            #最大値をとるときのルートを記憶する
            if scores[i][j] == scores[i-1][j]-1:
                how[i][j] = 1
            elif scores[i][j] == scores[i][j-1]-1:
                how[i][j] = 2
            else:
                how[i][j] = 3


    #最大値をとるときのルートを元にアラインメントを求める
    #ゴールから順にアラインメントし、全て終わったらループを抜ける
    col = len(main_seq)
    row = len(sub_seq)
    while(1):
        if how[col][row] == 1:
            changed_sub.append('-')
            changed_main.append(main_seq[col-1])
            col -= 1
        elif how[col][row] == 2:
            changed_main.append('-')
            changed_sub.append(sub_seq[row-1])
            row -= 1
        elif how[col][row] == 3:
            changed_sub.append(sub_seq[row-1])
            changed_main.append(main_seq[col-1])
            col -= 1
            row -= 1
        elif col == 0 and row != 0:
            changed_main.append('-')
            changed_sub.append(sub_seq[row-1])
            row -= 1
        elif row == 0 and col != 0:
            changed_sub.append('-')
            changed_main.append(main_seq[col-1])
            col -= 1

        if col == 0 and row == 0:
            break

    #終わりから順にアラインメント後の配列を作ったため、
    #逆順にする
    changed_sub.reverse()
    changed_main.reverse()
    
    #アラインメント後の配列を返り値とする
    return changed_sub, changed_main


#アラインメントする本数を読み込む
number = int(input())
#与えられた遺伝子配列を記憶する配列を用意
pre_sequences = []

#アラインメントする遺伝子配列の読み込み
for i in range(number):
    aline = list(input())
    pre_sequences.append(aline)

#ペア間の類似性スコアを計算する

#ペア間の類似性スコアを保存する配列を作成
similarity_scores = [[0 for j in range(number)] 
                for i in range(number)]


#全てのペアの類似性スコアを計算し、配列に保存する
for i in range(0,number-1):
    for j in range(i+1,number):
        similarity_scores[i][j] = calculate(pre_sequences[i], pre_sequences[j])
        similarity_scores[j][i] = calculate(pre_sequences[i], pre_sequences[j])


#各遺伝子配列の他の配列に対する類似性スコアの総和を求める
scores_sum = [0 for i in range(number)]

for i in range(number):
    scores_sum[i] = sum(similarity_scores[i])

#類似性スコアの一番大きい遺伝子配列を求める
#まず最も類似性スコアの大きい配列を1番目の配列とする
max_seq = 1
max_score = scores_sum[0]

#もし1番目より大きい類似性スコアをもつ配列があれば
#最大類似性スコアとその配列番号を更新する
for i in range(1,number):
    if scores_sum[i] >= max_score:
        max_score = scores_sum[i]
        max_seq = i+1

#最大類似性スコアをもつ配列を中心として,アラインメントを融合していく
ans_sequences = []
tmp_save = []
bar_save = []

for i in range(number):
    ans_seq = alignment(pre_sequences[max_seq-1], pre_sequences[i])
    for j in range(len(ans_seq[1])):
        if ans_seq[1][j] == '-' and ans_seq[0][j] != '-':
            tmp_save.append(j) 
    bar_save.insert(i, tmp_save)
    tmp_save = []
    ans_sequences.append(ans_seq[0])

#アラインメントする際に挿入された'-'を他の配列にも適応
for i in range(number):
    for j in range(len(bar_save[i])):
        for k in range(number):
            if k != i:
                ans_sequences[k].insert(bar_save[i][j], '-')


#結果のスコアを計算する
answer = 0
for i in range(len(ans_sequences[0])):
    for j in range(number-1):
        for k in range(j+1, number):
            answer += pair_score(ans_sequences[j][i], ans_sequences[k][i])

#作成したアラインメントを表示
print("\nmultiple alignment")
for i in ans_sequences:
    print(*i)

print("score:", answer)