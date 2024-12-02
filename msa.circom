pragma circom 2.1.8;

include "comparators.circom";
include "gates.circom";

/*---------------------------------------------------------------------
  alignment and sequences
---------------------------------------------------------------------*/                              
template T1() {
  signal input e, r, c;  // enabled, row, col
  signal output es, ese, ee;  
  signal e1  <== IsEqual()([e,1]); 
  signal rc  <== IsEqual()([c,r]);   
  signal c0  <== IsEqual()([c,0]);  
  signal r0  <== IsEqual()([r,0]);  
  signal nr0 <== NOT()(r0);  

  signal t1 <== rc * nr0;
  signal t2 <== c0 * nr0;

  es  <== e1 * r0;  
  ese <== e1 * t1;
  ee  <== e1 * t2;
}

template T2() {
  signal input e, c;  // enabled, col
  signal output ee;   // en east
  signal e1 <== IsEqual()([e,1]); 
  signal c0 <== IsEqual()([c,0]);  
  ee <== e1 * c0;
}

template pair_aln_seq(n, m) {  // n is seq_len, m is aln_len
  signal input seq[n];
  signal input aln[m];
  m = m + 1;
  signal output out;  

  component a[n][m];
 
  // first column
  a[0][0] = T1();
  a[0][0].e <== 1;
  a[0][0].r <== seq[0];
  a[0][0].c <== aln[0];
  for (var i = 1; i < n; i++) {
    a[i][0] = T1();
    a[i][0].e <== 0;
    // does not matter
    a[i][0].r <== seq[i];
    a[i][0].c <== aln[0];
  }
  // other column
  for (var i = 0; i < n; i++) {
    for (var j = 1; j < m; j++) {
      a[i][j] = T1();
      a[i][j].e <== a[i][j-1].ee + (i==0?0:a[i-1][j-1].ese) + (i==0?0:a[i-1][j].es);
      a[i][j].r <== seq[i];
      a[i][j].c <== j == (m-1) ? 0 : aln[j];   // pad with -
    }
  }

  component b[m];
  for (var i = 0; i < m; i++) {
    b[i] = T2();
    b[i].e <== i == 0 ? 0 : b[i-1].ee + a[n-1][i].es + a[n-1][i-1].ese;
    b[i].c <== i == (m-1) ? 0 : aln[i];
  }
  
  out <== b[m-1].ee;
}
template check_aln_seq(nseq, seq_len, aln_len) {
  signal input seq[nseq][seq_len];
  signal input aln[nseq][aln_len];
  signal output y;
  signal t[nseq];
  signal b[nseq];

  for (var i = 0; i < nseq; i++) {
    t[i] <== pair_aln_seq(seq_len, aln_len)(seq[i], aln[i]);
  }
  for (var i = 0; i < nseq; i++) {
    b[i] <== AND()(i==0?1:b[i-1], t[i]);
  }
  y <== b[nseq - 1];
}

/*---------------------------------------------------------------------
  alignment & score
---------------------------------------------------------------------*/
template scoring_system() {
  signal input x[2];
  signal output y;

  signal xeq  <== IsEqual()(x);   // [x[0],x[1]]);
  signal xneq <== NOT()(xeq);
  signal gap  <== OR()(IsEqual()([x[0],0]), IsEqual()([x[1],0]));
  signal bgap <== AND()(xeq,gap);   // both are gap
  signal ngap <== NOT()(gap);         // no gaps
  signal eq_ngap <== AND()(xeq, ngap);

  //    mismatch       both gap      equal but not gap
  y <== (xneq * -1) + (bgap * -1) + (eq_ngap * 1);
}

template pair_score(aln_len) {
  signal input aln[2][aln_len];
  signal output score;
  component c[aln_len];
  signal s[aln_len];

  for (var i = 0; i < aln_len; i++) {
    c[i] = scoring_system();
    c[i].x[0] <== aln[0][i];
    c[i].x[1] <== aln[1][i];
    s[i] <== (i == 0 ? 0 : s[i-1]) + c[i].y;
  }
  score <== s[aln_len - 1];
}

template msa_score(nseq, aln_len) {
  signal input aln[nseq][aln_len];
  signal output score;

  component ps[nseq][nseq];
  signal t[nseq * nseq];
  var count = 0;
  t[count] <== 0;
  for (var i = 0; i < nseq; i++) {
    for (var j = i; j < nseq; j++) {
      ps[i][j] = pair_score(aln_len);
      ps[i][j].aln[0] <== aln[i];
      ps[i][j].aln[1] <== aln[j];
      count++;
      t[count] <== t[count-1] + ps[i][j].score;
    }
  }
  score <== t[count];
}

template check_aln_score(nseq, aln_len) {
  signal input aln[nseq][aln_len];
  signal input score;
  signal output y;
  signal s <== msa_score(nseq, aln_len)(aln);
  y <== IsEqual()([s, score]);
}
/*---------------------------------------------------------------------
  main
---------------------------------------------------------------------*/                         
template Main(nseq, seq_len, aln_len) {
  signal input seq[nseq][seq_len];
  signal input aln[nseq][aln_len];
  signal input score;
  signal output y;

  signal aln_seq   <== check_aln_seq(nseq, seq_len, aln_len)(seq, aln);
  //y <== aln_seq;
  signal aln_score <== check_aln_score(nseq, aln_len)(aln, score);
  y <== AND()(aln_score, aln_seq);
}

component main {public [seq, score]} = Main(5,10,20);
