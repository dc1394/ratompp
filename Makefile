#
# プログラム名
#
PROG = ratompp

#
# ソースコードが存在する相対パス
#
VPATH = src/Alg src/ExCorr src/Fem1D src/Integ src/Ks src/Util

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
CXX = clang++

#
# C++コンパイラに与える、（最適化等の）オプション
#
CXXFLAGS = -Wall -Wextra -std=c++17 -O3 -m64 -I${MKLROOT}/include 

#
# リンク対象に含めるライブラリの指定
#
LDFLAGS = -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_intel_thread -lmkl_core \
		  -liomp5 -lpthread -lm -ldl -lxc
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
		$(CXX) $(LDFLAGS) $(CXXFLAGS) -o $@ $^

#
# プログラムのコンパイル
#
%.o: %.cpp
		$(CXX) $(CXXFLAGS) -c -MMD -MP $<

#
# make cleanの動作
#
clean:
		rm -f $(PROG) $(OBJS) $(DEPS)
