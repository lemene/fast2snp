# 1. 进入 src 目录编译 C++ 程序
# 2. 拷贝 scripts 中的 Python 脚本到 build/bin

.PHONY: all clean src scripts

# 默认目标
all: src scripts

# 编译 C++ 程序
src:
	$(MAKE) -C src

# 拷贝 Python 脚本
scripts:
	@mkdir -p build/bin
	@if [ -d scripts ]; then \
		find scripts -maxdepth 1 -name "*.py" -exec cp {} build/bin/ \; ; \
		echo "Python scripts copied to build/bin/"; \
	else \
		echo "Warning: scripts directory not found"; \
	fi

clean:
	$(MAKE) -C src clean
	rm -rf build

.PHONY: bin
bin: all