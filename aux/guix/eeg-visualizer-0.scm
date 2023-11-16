;; NOTE: This is untested as of yet.

(include "lsl.scm")

(define eeg-visualizer-0
  (package
    (name "eeg-visualizer-0")
    (version "0.0.1")
    (source (local-file "./../.." "eeg-visualizer-0" #:recursive #t))
    (build-system cmake-build-system)
    (arguments '(#:tests? #f))
    (inputs (list lsl sfml fftw))
    (home-page "https://github.com/publik-void/eeg-visualizer-0")
    (synopsis "Real-time EEG visualizer")
    (description #f)
    (license #f)))

eeg-visualizer-0
