%   本ファイルを, 「startup.m」という名前をつけて，カレントの起動ディレクトリに置く.
%   デフォルト設定では, $matlabroot\work. 
%   ここで $matlabrootは, MATLABファイルがインストールされたディレクトリである.
% 	
%   * なお, 次のコマンドにて, 初期設定を確認できる
%   >> get(0,'Factory')
%

%% figureウィンドウの 用紙設定：：A4サイズを基本とする
set(groot,'DefaultFigurePaperType','a4');
set(groot,'DefaultFigurePaperPosition',[4.15 10.08 12.69 9.52]);

%% figureウィンドウの 軸設定
set(groot,'defaultAxesBox','on');
set(groot,'defaultAxesLineWidth',1.5);		% 軸線の太さを1.5 ptにする
set(groot,'defaultAxesFontSize',28);		% フォントサイズを28 ptにする

set(groot,'defaultAxesFontName','Times New Roman');   % 軸ラベルのフォント種類を設定
% set(0,'defaultAxesFontName','Helvetica'); 	% Times New Romanの代わりに， HelveticaやTimesなども設定できる
% set(groot,'defaultAxesFontName','Times');

%% figureウィンドウの ライン設定
set(groot,'defaultLineLineWidth',2.0);	

%% figureウィンドウの フォント設定
set(groot,'defaultTextFontSize',28);	% グラフ内のフォントサイズを28 ptにする
set(groot,'defaultTextFontName','Times New Roman');
% set(0,'defaultTextFontName','Helvetica');
% set(groot,'defaultTextFontName','Times');

