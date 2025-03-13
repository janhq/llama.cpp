# Makefile for Cortex llamacpp engine - Build, Lint, Test, and Clean

CMAKE_EXTRA_FLAGS ?= ""
RUN_TESTS ?= false
LLM_MODEL_URL ?= "https://delta.jan.ai/tinyllama-1.1b-chat-v0.3.Q2_K.gguf"
EMBEDDING_MODEL_URL ?= "https://catalog.jan.ai/dist/models/embeds/nomic-embed-text-v1.5.f16.gguf"
CODE_SIGN ?= false
AZURE_KEY_VAULT_URI ?= xxxx
AZURE_CLIENT_ID ?= xxxx
AZURE_TENANT_ID ?= xxxx
AZURE_CLIENT_SECRET ?= xxxx
AZURE_CERT_NAME ?= xxxx
DEVELOPER_ID ?= xxxx

# Default target, does nothing
all:
	@echo "Specify a target to run"

# Build the Cortex engine
build-lib:
ifeq ($(OS),Windows_NT)
	@powershell -Command "cmake -B build $(CMAKE_EXTRA_FLAGS) -DLLAMA_BUILD_TESTS=OFF;"
	@powershell -Command "cmake --build build --config Release -j4 --target llama-server;"
else ifeq ($(shell uname -s),Linux)
	@cmake -B build $(CMAKE_EXTRA_FLAGS) -DLLAMA_BUILD_TESTS=OFF
	@cmake --build build --config Release -j4 --target llama-server;
else
	@cmake -B build $(CMAKE_EXTRA_FLAGS) -DLLAMA_BUILD_TESTS=OFF
	@cmake --build build --config Release -j4 --target llama-server;
endif

pre-package:
ifeq ($(OS),Windows_NT)
	@powershell -Command "mkdir -p build\bin; cp build\bin\llama-server.exe build\bin\;"
else ifeq ($(shell uname -s),Linux)
	@mkdir -p build/bin; \
	cp build/bin/llama-server build/bin/;
else
	@mkdir -p build/bin; \
	cp build/bin/llama-server build/bin/;
endif

codesign:
ifeq ($(CODE_SIGN),false)
	@echo "Skipping Code Sign"
	@exit 0
endif

ifeq ($(OS),Windows_NT)
	@powershell -Command "dotnet tool install --global AzureSignTool;"
	@powershell -Command 'azuresigntool.exe sign -kvu "$(AZURE_KEY_VAULT_URI)" -kvi "$(AZURE_CLIENT_ID)" -kvt "$(AZURE_TENANT_ID)" -kvs "$(AZURE_CLIENT_SECRET)" -kvc "$(AZURE_CERT_NAME)" -tr http://timestamp.globalsign.com/tsa/r6advanced1 -v ".\build\bin\llama-server.exe";'
else ifeq ($(shell uname -s),Linux)
	@echo "Skipping Code Sign for linux"
	@exit 0
else
	find "llama" -type f -exec codesign --force -s "$(DEVELOPER_ID)" --options=runtime {} \;
endif

package:
ifeq ($(OS),Windows_NT)
	@powershell -Command "7z a -ttar temp.tar build\bin\*; 7z a -tgzip llama.tar.gz temp.tar;"
else ifeq ($(shell uname -s),Linux)
	@tar -czvf llama.tar.gz build\bin;
else
	@tar -czvf llama.tar.gz build\bin;
endif

# run-e2e-test:
# ifeq ($(RUN_TESTS),false)
# 	@echo "Skipping tests"
# 	@exit 0
# endif
# ifeq ($(OS),Windows_NT)
# 	@powershell -Command "mkdir -p examples\server\build\engines\llama; cd examples\server\build; cp ..\..\..\build\engine.dll engines\llama; ..\..\..\.github\scripts\e2e-test-server-windows.bat server.exe $(LLM_MODEL_URL) $(EMBEDDING_MODEL_URL);"
# else ifeq ($(shell uname -s),Linux)
# 	@mkdir -p examples/server/build/engines/llama; \
# 	cd examples/server/build/; \
# 	cp ../../../build/libengine.so engines/llama/; \
# 	chmod +x ../../../.github/scripts/e2e-test-server-linux-and-mac.sh && ../../../.github/scripts/e2e-test-server-linux-and-mac.sh ./server $(LLM_MODEL_URL) $(EMBEDDING_MODEL_URL);
# else
# 	@mkdir -p examples/server/build/engines/llama; \
# 	cd examples/server/build/; \
# 	cp ../../../build/libengine.dylib engines/llama/; \
# 	chmod +x ../../../.github/scripts/e2e-test-server-linux-and-mac.sh && ../../../.github/scripts/e2e-test-server-linux-and-mac.sh ./server $(LLM_MODEL_URL) $(EMBEDDING_MODEL_URL);
# endif

# run-e2e-submodule-test:
# ifeq ($(RUN_TESTS),false)
# 	@echo "Skipping tests"
# 	@exit 0
# endif
# ifeq ($(OS),Windows_NT)
# 	@powershell -Command "python -m pip install --upgrade pip"
# 	@powershell -Command "python -m pip install requests;"
# 	@powershell -Command "mkdir -p examples\server\build\engines\llama; cd examples\server\build; cp ..\..\..\build\engine.dll engines\llama; python ..\..\..\.github\scripts\e2e-test-server.py server $(LLM_MODEL_URL) $(EMBEDDING_MODEL_URL);"
# else ifeq ($(shell uname -s),Linux)
# 	python -m pip install --upgrade pip;
# 	python -m pip install requests;
# 	@mkdir -p examples/server/build/engines/llama; \
# 	cd examples/server/build/; \
# 	cp ../../../build/libengine.so engines/llama/; \
# 	python  ../../../.github/scripts/e2e-test-server.py server $(LLM_MODEL_URL) $(EMBEDDING_MODEL_URL);
# else
# 	python -m pip install --upgrade pip;
# 	python -m pip install requests;
# 	@mkdir -p examples/server/build/engines/llama; \
# 	cd examples/server/build/; \
# 	cp ../../../build/libengine.dylib engines/llama/; \
# 	python  ../../../.github/scripts/e2e-test-server.py server $(LLM_MODEL_URL) $(EMBEDDING_MODEL_URL);
# endif

# run-e2e-weekend-test:
# ifeq ($(RUN_TESTS),false)
# 	@echo "Skipping tests"
# 	@exit 0
# endif
# ifeq ($(OS),Windows_NT)
# 	@powershell -Command "python -m pip install --upgrade pip"
# 	@powershell -Command "python -m pip install requests;"
# 	@powershell -Command "mkdir -p examples\server\build\engines\llama; cd examples\server\build; cp ..\..\..\build\engine.dll engines\llama; python ..\..\..\.github\scripts\e2e-test-server-weekend.py server;"
# else ifeq ($(shell uname -s),Linux)
# 	python -m pip install --upgrade pip;
# 	python -m pip install requests;
# 	@mkdir -p examples/server/build/engines/llama; \
# 	cd examples/server/build/; \
# 	cp ../../../build/libengine.so engines/llama/; \
# 	python  ../../../.github/scripts/e2e-test-server-weekend.py server;
# else
# 	python -m pip install --upgrade pip;
# 	python -m pip install requests;
# 	@mkdir -p examples/server/build/engines/llama; \
# 	cd examples/server/build/; \
# 	cp ../../../build/libengine.dylib engines/llama/; \
# 	python  ../../../.github/scripts/e2e-test-server-weekend.py server;
# endif