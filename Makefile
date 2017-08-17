SHELL := /nix/store/jjg5zsfjl35b239f6znp0nfa9775qy2a-bash-4.4-p12/bin/bash

print-mesh:
	cargo build --target=wasm32-unknown-emscripten --release --bin print-mesh
	mkdir -p site
	find target/wasm32-unknown-emscripten/release/deps -type f -name "*.wasm" | xargs -I {} cp {} site/site.wasm
	find target/wasm32-unknown-emscripten/release/deps -type f ! -name "*.asm.js" -name "*.js" | xargs -I {} cp {} site/site.js