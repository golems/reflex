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


(defparameter *test-frames*
  `((frame :name frame_0
           :configuration q0)
    (frame :name frame_1
           :configuration q1)
    (frame :name frame_2
           :configuration q2)))

(defun prefix-frame-1 (prefix-parent prefix frame)
  (cons :frame
        (loop for k = (cdr frame) then (cddr k)
           while k
           for key = (first k)
           for val = (second k)
           append (list key
                        (cond
                          ((or (eq :name key) (eq :configuration key)
                               (and (eq :parent key)
                                    (not (equal val prefix-parent))))
                           (concatenate 'string prefix (string val)))
                          (t
                           val))))))

(defun prefix-frames (prefix-parent prefix frames)
  (loop for f in frames
       collect (prefix-frame-1 prefix-parent prefix f)))


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
  (emit-indices stream
                (loop
                   with i = -1
                   for frame in frames
                   for config = (frame-configuration frame)
                   when config
                   collect (list config (incf i)))
                max))

(defun float-string (f)
  (map 'string (lambda (c)
                 (if (or (eq c #\d)
                         (eq c #\D))
                     #\E
                     c))
                 (format nil "~F" f)))

(defun emit-fixed-frame (frame ptr stream)
  (let ((q (frame-quaternion frame))
        (v (frame-translation frame)))
    (dotimes (i 4)
      (format stream "~&    (~A)[~D] = ~A;"
              ptr i (float-string (elt q i))))
    (dotimes (i 3)
      (format stream "~&    (~A)[~D] = ~A;"
              ptr (+ 4 i) (float-string (elt v i))))))

(defun emit-revolute-frame (frame ptr stream)
  (let ((v (frame-translation frame))
        (axis (frame-axis frame)))
    (format stream "~&    {")
    (format stream "~&        static const double axis[3] = {~A,~A,~A};"
            (float-string (elt axis 0))
            (float-string (elt axis 1))
            (float-string (elt axis 2)))
    (format stream "~&        aa_tf_axang2quat( axis, q[~A], (~A) );"
            (frame-configuration frame) ptr )
    (format stream "~&    }")
    (dotimes (i 3)
      (format stream "~&    (~A)[~D] = ~A;"
              ptr (+ 4 i) (float-string (elt v i))))))

(defun emit-parents-array (name frames &key (stream t))
  (format stream "ssize_t ~A[] = {~&~{~&    ~A~^,~}~&};"
          name
          (loop for f in frames
             for p = (frame-parent f)
             collect (if p p -1))))

(defun emit-rel-fun (name frames &key (stream t))
  (format stream "void ~A( double *q, double *e ) ~&{~&" name)
  (loop
     for frame in frames
     for i from 0
     for ptr = (format nil "e+7*~D" i)
     do
       (format stream "~&    // Frame: ~A, type: ~A"
               (frame-name frame)
               (frame-type frame))
       (case (frame-type frame)
         (:revolute (emit-revolute-frame frame ptr stream))
         (:fixed (emit-fixed-frame frame ptr stream))
         (otherwise (error "Unknown type of frame ~A" (frame-name frame)))))
  (format stream "~&}"))

(defun write-frame-header (filename frames &key
                           frame-max
                           configuration-max)
  (with-open-file (f filename :direction :output :if-exists :supersede)
    (format f "~&/* AUTOGENERATED BY REFLEX FRAME GENERATOR */")
    (format f "~&~&/* NAME INDICES */")
    (emit-name-indices frames :stream f)
    (format f "~&~&/* CONFIGURATION INDICES */")
    (emit-config-indices frames :stream f)
    ) )

;; (defun write-frame-source (name frames)
;;   (with-open-file (f name :direction :output :if-exists :supersede)
;;     (format f "~&/* AUTOGENERATED BY REFLEX FRAME GENERATOR */")
