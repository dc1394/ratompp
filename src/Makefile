#
# プログラム名
#
PROG = ratompp

#
# プログラムをインストールする場所の指定
#
DESTDIR = ../

#
# ソースコードが存在する相対パス
#
VPATH = Alg ExCorr Fem1D Integ Ks Util

#
# コンパイル対象のソースファイル群（カレントディレクトリ以下の*.cppファイル）
#
SRCS = $(shell find * -name "*.cpp")

#
# ターゲットファイルを生成するために利用するオブジェクトファイル
#
OBJDIR = 
ifeq "$(strip $(OBJDIR))" ""
  OBJDIR = .
endif

OBJS = $(addprefix $(OBJDIR)/, $(notdir $(SRCS:.cpp=.o)))

#
# *.cppファイルの依存関係が書かれた*.dファイル
#
DEPS = $(OBJS:.o=.d)

#
# C++コンパイラの指定
#
CXX = g++

#
# C++コンパイラに与える、（最適化等の）オプション
#
CXXFLAGS = -Wall -Wextra -std=c++17 -O3 -mtune=native -march=native

#
# リンク対象に含めるライブラリの指定
#
LDFLAGS = -lopenblas -lpthread -lm -ldl -lxc
#
# makeの動作
#
all: $(PROG) ;

#
# 依存関係を解決するためのinclude文
#
-include $(DEPS)

#
# プログラムのリンク
#
$(PROG): $(OBJS)
		$(CXX) $^ $(LDFLAGS) $(CXXFLAGS) -o $@

#
# プログラムのコンパイル
#
%.o: %.cpp
		$(CXX) $(CXXFLAGS) -c -MMD -MP $<


#
# プログラムのインストール
#
install: $(PROG)
		strip $(PROG)
		cp $(PROG) $(DESTDIR)$(PROG)

#
# make cleanの動作
#
clean:
		rm -f $(PROG) $(OBJS) $(DEPS)
