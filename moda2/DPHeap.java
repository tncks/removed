package moda;

import java.util.LinkedList;

import msutil.PGraph;
import msutil.Scoring;

public class DPHeap extends LinkedList<DPPeptide> {
    // TODO: check suggestion at the bottom line

    static final int Capacity = 100; // max heap capacity //MODplus

    public DPHeap() {
        for (int i = 0; i < Capacity; i++)
            this.add(new DPPeptide());
    }

    public void setStemNo(int stem) {
        for (DPPeptide dp : this)
            dp.stem = stem;

        //for (int i = 0; i < this.size(); i++) this.get(i).stem = stem;
    }

    public boolean isConfident() {
        for (int i = 0; i < this.size() && i < 10; i++) {
            if (this.get(i).isConfident()) return true;
        }
        return false;
    }

    public boolean insert(DPPeptide dp){

        if( dp.score < 1 || dp.compareTo( this.getLast() ) == 1 ) return false;

        int i=0;
        for( DPPeptide x : this ) {
            if( dp.isSame(x) ) break;
            if( x.compareTo(dp) == 1 ){
                this.add(i, dp);
                this.removeLast();
                break;
            }
            i++;
        }
        return true;
    }


    public boolean insertInternal(DPPeptide dp) {
        if (dp.score < 1 || dp.compareTo(this.getLast()) == 1) return false;

        int i = 0;
        for (DPPeptide x : this) {
            if (dp.isSame(x)) break;
            if (x.compareTo(dp) == 1) {
                this.add(i, dp);
                this.removeLast();
                break;
            }
            i++;
        }
        return true;
    }

    public void insertAll(DPHeap heap) {
        for (DPPeptide x : heap) {
            insertInternal(x);
        }

    }


    public int evaluate(PGraph graph) { //for moda
        int i;
        int maxScore = 0;
        DPPeptide ele;
        for (i = 0; i < this.size(); i++) {
            ele = this.get(i);
            if (ele.score < 1) {
                this.remove(i);
                i--;
                continue;
            }
            this.get(i).evaluatePSM(graph);
            if (maxScore < this.get(i).score) maxScore = this.get(i).score;
        }

        for (i = 0; i < this.size(); i++) {
            if (this.get(i).score == maxScore) {
                this.get(i).score -= 1;
            }
        }
        this.sort(new DPPeptideRefinedComparator());
        if (this.size() == 0 || this.get(0).score < 1) return 0;
        return i;
    }

    public int reArr(PGraph graph) throws NullPointerException { //for modplus
        int i;
        int maxScore = 0;
        DPPeptide ele;
        for (i = 0; i < this.size(); i++) {
            ele=this.get(i);
            if (ele.score < 1) {
                this.remove(i);
                i--;
                continue;
            }
            this.get(i).score = Scoring.getModARankScore(ele.peptide, ele.ptms, graph);
            ele=this.get(i);
            if (maxScore < ele.score) maxScore = ele.score;
        }

        this.sort(new DPPeptideRefinedComparator());
        int cut = maxScore / 2;
        for (i = 0; i < this.size(); i++) {
            if (this.get(i).score < cut) this.removeRange(i, this.size());
        }

        if (this.size() == 0 || this.get(0).score < 1) return 0;
        return i;
    }

}

/**
 *
 *
 *
 */
// Suggestions (advice for code clean up) - writing time based on this date: 2025 03 18 dawn (12:15 am)

/*
current

G1

G2



public int reArr(PGraph graph) throws NullPointerException { //for modplus
        int i = 0;
        int maxScore = 0;
        for (i = 0; i < this.size(); i++) {
            if (this.get(i).score < 1) {
                this.remove(i);
                i--;
                continue;
            }
            this.get(i).score = Scoring.getModARankScore(this.get(i).peptide, this.get(i).ptms, graph);
            if (maxScore < this.get(i).score) maxScore = this.get(i).score;
        }

        this.sort(new DPPeptideRefinedComparator());
        int cut = maxScore / 2;
        for (i = 0; i < this.size(); i++) {
            if (this.get(i).score < cut) this.removeRange(i, this.size());
        }

        if (this.size() == 0 || this.get(0).score < 1) return 0;
        return i;
    }


public int reArr(PGraph graph) throws NullPointerException {
    int maxScore = this.stream()
            .filter(dp -> dp.score >= 1)  // score가 1 이상인 것만 필터링
            .peek(dp -> dp.score = Scoring.getModARankScore(dp.peptide, dp.ptms, graph))  // 점수 재계산
            .mapToInt(dp -> dp.score)  // 점수만 추출하여
            .max()  // 최대 점수 찾기
            .orElse(0);  // 빈 리스트의 경우 0 반환

    // score가 최대 점수와 같은 peptide들의 점수에서 1을 차감
    this.stream()
            .filter(dp -> dp.score == maxScore)
            .forEach(dp -> dp.score -= 1);

    // 정렬 후 일부 제거
    this.sort(new DPPeptideRefinedComparator());
    this.removeIf(dp -> dp.score < maxScore / 2);  // 최대 점수의 절반 미만은 제거

    // 최종적으로 남은 peptide 수를 반환
    return (this.size() == 0 || this.get(0).score < 1) ? 0 : this.size();
}

public int reArr(PGraph graph) throws NullPointerException {
    if (isEmpty()) return 0;

    List<DPPeptide> updatedList = stream()
            .filter(dp -> dp.score >= 1)
            .peek(dp -> dp.score = Scoring.getModARankScore(dp.peptide, dp.ptms, graph))
            .sorted(new DPPeptideRefinedComparator())
            .collect(Collectors.toList());

    if (updatedList.isEmpty()) return 0;

    int maxScore = updatedList.get(0).score;
    int cut = maxScore / 2;

    clear();
    addAll(updatedList.stream()
            .filter(dp -> dp.score >= cut)
            .collect(Collectors.toList()));

    return size();
}
 */
