.PHONY: all build install clean test deps

# Binary name
BINARY=riskdays_go
INSTALL_PATH=/usr/local/bin

# Build the binary
all: build

# Download dependencies
deps:
	go mod download
	go mod tidy

# Build the binary
build: deps
	mkdir -p bin
	go build -o bin/$(BINARY) main.go

# Install binary to system path
install: build
	cp bin/$(BINARY) $(INSTALL_PATH)/$(BINARY)
	chmod +x $(INSTALL_PATH)/$(BINARY)

# Run tests
test:
	go test -v ./...

# Run tests with coverage
coverage:
	go test -v -coverprofile=coverage.out ./...
	go tool cover -html=coverage.out

# Clean build artifacts
clean:
	rm -rf bin/
	rm -f coverage.out

# Build for multiple platforms
build-all: deps
	mkdir -p bin
	GOOS=darwin GOARCH=amd64 go build -o bin/$(BINARY)-darwin-amd64 main.go
	GOOS=darwin GOARCH=arm64 go build -o bin/$(BINARY)-darwin-arm64 main.go
	GOOS=linux GOARCH=amd64 go build -o bin/$(BINARY)-linux-amd64 main.go
	GOOS=windows GOARCH=amd64 go build -o bin/$(BINARY)-windows-amd64.exe main.go

# Run the binary
run: build
	./bin/$(BINARY)
