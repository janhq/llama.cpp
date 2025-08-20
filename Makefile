<<<<<<< HEAD
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
NOTARIZE ?= false
QUILL_NOTARY_KEY_ID ?= ""
QUILL_NOTARY_ISSUER ?= ""
QUILL_NOTARY_KEY ?= "/tmp/notary-key.p8"

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
	find "build/bin" -type f -exec codesign --force -s "$(DEVELOPER_ID)" --options=runtime {} \;
endif

notarize:
ifeq ($(NOTARIZE),false)
	@echo "Skipping Notarization"
	@exit 0
endif

ifeq ($(OS),Windows_NT)
	@echo "Skipping Notarization for Windows"
	@exit 0
else ifeq ($(shell uname -s),Linux)
	@echo "Skipping Notarization for Linux"
	@exit 0
else
	@echo "Starting notarization for macOS binaries..."
	@find build/bin -type f -perm +111 -exec env QUILL_NOTARY_KEY_ID="$(QUILL_NOTARY_KEY_ID)" QUILL_NOTARY_ISSUER="$(QUILL_NOTARY_ISSUER)" QUILL_NOTARY_KEY="$(QUILL_NOTARY_KEY)" quill notarize {} \;
	@echo "All macOS binaries notarized successfully"
endif

package:
ifeq ($(OS),Windows_NT)
	@powershell -Command "7z a -ttar temp.tar build\bin\*; 7z a -tgzip llama.tar.gz temp.tar;"
else ifeq ($(shell uname -s),Linux)
	@tar -czvf llama.tar.gz build/bin;
else
	@tar -czvf llama.tar.gz build/bin;
endif
=======
define newline


endef

$(error Build system changed:$(newline)\
The Makefile build has been replaced by CMake.$(newline)$(newline)\
For build instructions see:$(newline)\
https://github.com/ggml-org/llama.cpp/blob/master/docs/build.md$(newline)${newline})
>>>>>>> 37f10f95 (make : remove make in favor of CMake (#15449))
