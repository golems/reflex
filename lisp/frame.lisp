;;;; -*- mode: lisp -*-
;;;;
;;;; Copyright (c) 2013, Georgia Tech Research Corporation
;;;; All rights reserved.
;;;;
;;;; Author(s): Neil T. Dantam <ntd@gatech.edu>
;;;; Georgia Tech Humanoid Robotics Lab
;;;; Under Direction of Prof. Mike Stilman <mstilman@cc.gatech.edu>
;;;;
;;;;
;;;; This file is provided under the following "BSD-style" License:
;;;;
;;;;
;;;;   Redistribution and use in source and binary forms, with or
;;;;   without modification, are permitted provided that the following
;;;;   conditions are met:
;;;;
;;;;   * Redistributions of source code must retain the above copyright
;;;;     notice, this list of conditions and the following disclaimer.
;;;;
;;;;   * Redistributions in binary form must reproduce the above
;;;;     copyright notice, this list of conditions and the following
;;;;     disclaimer in the documentation and/or other materials provided
;;;;     with the distribution.
;;;;
;;;;   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
;;;;   CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
;;;;   INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
;;;;   MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
;;;;   DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
;;;;   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
;;;;   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
;;;;   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
;;;;   USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
;;;;   AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
;;;;   LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
;;;;   ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
;;;;   POSSIBILITY OF SUCH DAMAGE.

(in-package :reflex)

(defmacro def-frame-var (var doc)
  (let ((frame (gensym "FRAME")))
    `(defun ,(intern (concatenate 'string "FRAME-" (string var)))
         (,frame)
       ,doc
       (destructuring-bind (,frame &key ,var &allow-other-keys) ,frame
         (assert (eq :frame ,frame))
         ,var))))

(def-frame-var name "name of the frame")
(def-frame-var type "type of the frame, fixed or revolute")
(def-frame-var parent "parent of this frame")
(def-frame-var axis "frame axis of rotation, (x y z)")
(def-frame-var configuration "name of frame configuration")
(def-frame-var quaternion "fixed transform quaternion")
(def-frame-var translation "translation (x y z)")
(def-frame-var offset "Configuration offset")


(defparameter *test-frames*
  `((frame :name frame_0
           :configuration q0)
    (frame :name frame_1
           :configuration q1)
    (frame :name frame_2
           :configuration q2)))

(defun prefix-frame-1 (prefix frame parent)
  (labels ((prefix (val) (concatenate 'string prefix (string val))))
    (append (list :frame
                  :parent (if (and (frame-parent frame)
                                   (not (equal parent (frame-parent frame))))
                              (prefix (frame-parent frame))
                              parent))
            (loop for k = (cdr frame) then (cddr k)
               while k
               for key = (first k)
               for val = (second k)
               unless (eq key :parent)
               append (list key
                            (cond
                              ((or (eq :name key) (eq :configuration key))
                               (prefix val))
                              (t
                               val)))))))

(defun prefix-frames (parent prefix frames)
  (loop for f in frames
       collect (prefix-frame-1 prefix f parent)))


(defun make-fixed-frame (name parent orientation
                         &key
                         (x 0)
                         (y 0)
                         (z 0))
  `(:frame :name ,name
          :parent ,parent
          :type :fixed
          :translation (,x ,y ,z)
          :quaternion ,(let ((q (amino:quaternion orientation)))
                            (list (amino::quaternion-x q)
                                  (amino::quaternion-y q)
                                  (amino::quaternion-z q)
                                  (amino::quaternion-w q)))))

(defun emit-indices (stream indices max)
  (format stream "~{~&#define ~{~A ~D~}~}"
          (if max
              (cons (list max (length indices)) indices)
              indices)))

(defun emit-name-indices (frames &key (stream t) max)
  (emit-indices stream
                (loop
                   for frame in frames
                   for i from 0
                   collect (list (frame-name frame)
                                 i))
                max))


(defun emit-config-indices (frames &key (stream t) max)
  (let ((vars (remove-duplicates (loop
                                    for frame in frames
                                    for config = (frame-configuration frame)
                                    when config
                                    collect config)
                                 :test #'equal)))
    (emit-indices stream
                  (loop
                     for i from 0
                     for v in vars
                     collect (list v i))
                  max)))

(defun float-string (f)
  (if (numberp f)
      (map 'string (lambda (c)
                     (if (or (eq c #\d)
                             (eq c #\D))
                         #\E
                         c))
           (format nil "~F" f))
      f))

(defun emit-fixed-frame (frame stream e)
  (let ((q (frame-quaternion frame))
        (v (frame-translation frame))
        (ptr (format nil "(~A+7*~A)" e (frame-name frame))))
    (dotimes (i 4)
      (format stream "~&    ~A[~D] = ~A;"
              ptr i (float-string (elt q i))))
    (dotimes (i 3)
      (format stream "~&    ~A[~D] = ~A;"
              ptr (+ 4 i) (float-string (elt v i))))))

(defun emit-quat (stream ptr angle x y z)
  (let* ((n (sqrt (+ (* x x) (* y y) (* z z))))
         (x (/ x n))
         (y (/ y n))
         (z (/ z n)))
    (cond
      ((= x 1)  (format stream "~&    aa_tf_xangle2quat( ~A, ~A );" angle ptr))
      ((= x -1) (format stream "~&    aa_tf_xangle2quat( -~A, ~A );" angle ptr))
      ((= y 1)  (format stream "~&    aa_tf_yangle2quat( ~A, ~A );" angle ptr))
      ((= y -1) (format stream "~&    aa_tf_yangle2quat( -~A, ~A );" angle ptr))
      ((= z 1)  (format stream "~&    aa_tf_zangle2quat( ~A, ~A );" angle ptr))
      ((= z -1) (format stream "~&    aa_tf_zangle2quat( -~A, ~A );" angle ptr))
       (t
        (format stream "~&    {")
        (format stream "~&        static const double axis[3] = {~A,~A,~A};"
                (float-string x)
             (float-string y)
             (float-string z))
     (format stream "~&        aa_tf_axang2quat2( axis, ~A, ~A );"
             angle
             ptr)
     (format stream "~&    }")))))

(defun emit-revolute-frame (frame stream q-array e-array)
  (let ((v (frame-translation frame))
        (axis (frame-axis frame))
        (ptr (format nil "(~A+7*~A)" e-array (frame-name frame)))
        (q (format nil "~A[~A]" q-array (frame-configuration frame)))
        (offset (frame-offset frame)))
    (emit-quat stream ptr (if offset
                              (format nil "(~A+~A)" q offset)
                              q)
               (elt axis 0) (elt axis 1) (elt axis 2))
    (dotimes (i 3)
      (format stream "~&    ~A[~D] = ~A;"
              ptr (+ 4 i) (float-string (elt v i))))))

(defun emit-parents-array (name frames &key (stream t))
  (format stream "~&const ssize_t ~A[] = {~&~{~&    ~A~^,~}~&};"
          name
          (loop for f in frames
             for p = (frame-parent f)
             collect (if p p -1))))


(defun emit-names-array (name frames &key (stream t))
  (format stream "~&const char * ~A[] = {~&~{~&    \"~A\"~^,~}~&};"
          name
          (loop for f in frames
             for p = (frame-name f)
             collect (if p p -1))))

(defun emit-rel-fun (name frames &key (stream t))
  (format stream "~&void ~A~&( double *restrict q, double *restrict e ) ~&{~&" name)
  (loop
     for frame in frames
     for i from 0
     do
       (format stream "~&    // Frame: ~A, type: ~A"
               (frame-name frame)
               (frame-type frame))
       (case (frame-type frame)
         (:revolute (emit-revolute-frame frame stream "q" "e"))
         (:fixed (emit-fixed-frame frame stream "e"))
         (otherwise (error "Unknown type of frame ~A" (frame-name frame)))))
  (format stream "~&}"))


(defun emit-abs-fun (name frames &key normalize (stream t))
  (format stream "~&void ~A~&( const double * restrict rel, double * restrict abs ) ~&{~&" name)
  (loop
     for f in frames
     for name = (frame-name f)
     for rel = (format nil "rel + 7*~A" name)
     for abs = (format nil "abs + 7*~A" name)
     for parent-name = (frame-parent f)
     for par = (format nil "abs + 7*~A" parent-name)
     do
       (if parent-name
           (format stream "~&    aa_tf_qutr_~A( ~A, ~A, ~A );"
                   (if normalize "mulnorm" "mul")
                   par rel abs)
           (format stream "~&    memcpy( ~A, ~A, 7*sizeof(abs[0]) );"
                   abs rel)))
  (format stream "~&}"))

(defun emit-dot (frames &key (stream t))
  (format stream  "~&digraph {~&")
  (format stream  "~&   graph [rankdir=LR];")
  (loop for f in frames
     for name = (frame-name f)
     for parent = (frame-parent f)
     when parent
     do (format stream  "~&    ~A -> ~A;"
                parent name))
  (format stream  "~&}"))

(defun write-frame-files (header-file source-file frames &key
                          headers
                          system-headers
                          frame-max
                          configuration-max
                          dot-file
                          normalize
                          (relative-function "rel_tf")
                          (absolute-function "abs_tf")
                          (parents-array "parents")
                          (names-array "parents")
                          )
  ;; header
  (with-open-file (f header-file :direction :output :if-exists :supersede)
    (format f "~&/* AUTOGENERATED BY REFLEX FRAME GENERATOR */")
    (format f "~{~&#include <~A>~}" system-headers)
    (format f "~{~&#include \"~A\"~}~%" headers)
    (format f "~&extern const char *~A[];" names-array)
    (format f "~&extern const ssize_t ~A[];" parents-array)
    (format f "~&~%/* NAME INDICES */")
    (emit-name-indices frames :stream f :max frame-max)
    (format f "~&~%/* CONFIGURATION INDICES */")
    (emit-config-indices frames :stream f :max configuration-max)
    (format f "~&~%/* Compute Relative Transforms */")
    (format f "~&void ~A~&( double *AA_RESTRICT q, double *AA_RESTRICT e );"
            relative-function)
    (format f "~&~%/* Compute Absolute Transforms */")
    (format f "~&void ~A~&( const double * AA_RESTRICT rel, double * AA_RESTRICT abs );"
            absolute-function))
  ;; source
  (with-open-file (f source-file :direction :output :if-exists :supersede)
    (format f "~&/* AUTOGENERATED BY REFLEX FRAME GENERATOR */")
    (format f "~{~&#include <~A>~}" (list "amino.h" "reflex.h" header-file))
    (format f "~{~&#include \"~A\"~}" headers)
    (emit-rel-fun relative-function frames :stream f)
    (emit-parents-array parents-array frames :stream f)
    (emit-names-array names-array frames :stream f)
    (emit-abs-fun absolute-function frames :stream f :normalize normalize)
    )
  ;; dot
  (when dot-file
    (with-open-file (f dot-file :direction :output :if-exists :supersede)
      (emit-dot frames :stream f))))
