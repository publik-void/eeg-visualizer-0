(add-to-load-path (dirname (current-filename)))

(use-modules
  (guix packages)
  (guix gexp)
  (guix build-system cmake)
  (lsl)
  (gnu packages pkg-config)
  (gnu packages game-development)
  (gnu packages algebra))

;;; TODO: Add a font that the visualizer will use by default

(define eeg-visualizer-0
  (package
    (name "eeg-visualizer-0")
    (version "0.0.1")
    (source (local-file "./../.." "eeg-visualizer-0" #:recursive? #t))
    (build-system cmake-build-system)
    (arguments '(#:tests? #f))
    (native-inputs (list pkg-config))
    (inputs (list lsl sfml fftwf))
    (home-page "https://github.com/publik-void/eeg-visualizer-0")
    (synopsis "Real-time EEG visualizer")
    (description #f)
    (license #f)))

eeg-visualizer-0
